#!/usr/bin/env python
"""
A particle decay time simulator for investigation of the spectral add up method 

xaratustrah@github 2022-2025

"""

import argparse
from loguru import logger
import sys, os
import toml
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import curve_fit, OptimizeWarning
import numpy as np
from pydantic import BaseModel, Field
from tqdm import tqdm


# Change global font size for publication
plt.rcParams.update({"font.size": 14})

# handle scipy and numpy warnings as errors
warnings.simplefilter("error", OptimizeWarning)
np.seterr(all="ignore")  # raise


class Simulator:
    def __init__(self, config):
        # Extracting individual parameters
        self.params_tau_seed = config.params.tau_seed
        self.params_timestep = config.params.timestep
        self.params_n_sim_steps = config.params.n_sim_steps
        self.params_n_sim = config.params.n_sim
        self.params_n_trials = config.params.n_trials
        self.params_n_decay_steps = config.params.n_decay_steps
        self.params_mean_ion = config.params.mean_ion
        self.params_stdv_ion = config.params.stdv_ion
        self.params_mean_bkgnd = config.params.mean_bkgnd
        self.params_stdv_bkgnd = config.params.stdv_bkgnd
        self.params_empty_shots = config.params.empty_shots
        self.params_empty_shots_probability = config.params.empty_shots_probability

        # handle trailing slash and file path properly
        self.settings_output_path = os.path.join(config.settings.output_path, "")
        self.settings_plot_every_event = config.settings.plot_every_event
        self.settings_save_npz = config.settings.save_npz
        self.settings_tasks = config.settings.tasks
        self.settings_plot_titles = config.settings.plot_titles

        if not os.path.isdir(self.settings_output_path):
            logger.error("Output file path does not exist.")
            exit()

        logger.info("Simulation start. Enter command or ctrl-C to abort.")
        logger.info(f"Output file path: {self.settings_output_path}")

        self.simulation_error_flag = False

    @staticmethod
    def gaussian_function(x, *p):
        return p[0] * np.exp(-((x - p[1]) ** 2) / (2.0 * p[2] ** 2))

    @staticmethod
    def exponential_function(x, *p):
        # Exponential function + offset
        return p[0] * np.exp(-x / p[1]) + p[2]  # amplitude, tau and constant

    def create_from_events_boxcar(self, x):
        b_arr = np.array([])

        for i in range(self.params_n_sim):
            # create a single random number based on an exponential distribution
            num = np.random.exponential(self.params_tau_seed)

            idx = np.where(x < num)[0][-1]
            b = np.zeros(len(x))
            ii = np.r_[0:idx]
            b[ii] = self.params_mean_ion

            # decay happens in the same time bin
            if self.params_n_decay_steps == 0:
                jj = np.r_[idx : len(x)]
                b[jj] = self.params_mean_bkgnd

            # or latest in the next time bin
            elif self.params_n_decay_steps == 1:
                if not idx + 1 > len(x):
                    jj = np.r_[idx : idx + 1]
                    b[jj] = self.params_mean_ion / 2
                jj = np.r_[idx + 1 : len(x)]
                b[jj] = self.params_mean_bkgnd

            b_arr = np.append(b_arr, b)
            # plt.step(np.arange(len(b)), b, where = 'post')

        b_arr = np.reshape(b_arr, (self.params_n_sim, len(b)))
        b_arr_avg = np.average(b_arr, axis=0)
        return b_arr_avg

    def create_from_distribution(self, x):
        data = np.random.exponential(self.params_tau_seed, size=len(x))
        hist = np.histogram(
            data,
            bins=int(
                self.params_timestep * self.params_n_sim_steps / self.params_timestep
            ),

            # TODO:
            # the range of hist and x do not match with certain self.params_n_sim_steps
            range=(0, self.params_timestep * self.params_n_sim_steps),

        )[0]
        
        return hist / np.max(hist)

    def create_from_events_with_fluctuations(
        self, x, trial_number, add_empty_shots=False
    ):
        b_arr = np.array([])
        empty_shot_mask = np.random.choice(
            [0, 1],
            size=self.params_n_sim,
            p=[
                self.params_empty_shots_probability,
                1 - self.params_empty_shots_probability,
            ],
        )

        for i in range(self.params_n_sim):
            num = np.random.exponential(self.params_tau_seed)
            b = np.random.normal(self.params_mean_ion, self.params_stdv_ion, len(x))
            idx = np.where(x < num)[0][-1]

            if self.params_n_decay_steps == 0:
                b[idx : len(x)] = np.random.normal(
                    self.params_mean_bkgnd, self.params_stdv_bkgnd, len(x) - idx
                )

            elif self.params_n_decay_steps == 1:
                b[idx : idx + 1] = self.params_mean_ion / 2
                b[idx + 1 : len(x)] = np.random.normal(
                    self.params_mean_bkgnd, self.params_stdv_bkgnd, len(x) - idx - 1
                )

            if add_empty_shots and empty_shot_mask[i]:
                b = np.random.normal(
                    self.params_mean_bkgnd, self.params_stdv_bkgnd, len(x)
                )

            b_arr = np.append(b_arr, b)

            if self.settings_plot_every_event:
                x = np.arange(len(b)) * self.params_timestep
                y = b
                plt.step(x, y, where="post", color="#008080")  # color teal

                if self.settings_plot_titles:
                    plt.title(r'$\tau_{seed} =$'+str(self.params_tau_seed) + ' [s]')

                plt.xlabel("Time [s]")
                plt.ylabel("Amplitude [a.u.]")
                plt.grid()
                outfilename = f"{self.settings_output_path}trial{trial_number:04}_decay{i:04}_ts{self.params_tau_seed:.2e}"
                plt.tight_layout()
                plt.savefig(outfilename + ".png")
                plt.close()
                if self.settings_save_npz:
                    np.savez(outfilename + ".npz", x=x, y=y)

        b_arr = np.reshape(b_arr, (self.params_n_sim, len(b)))
        b_arr_avg = np.average(b_arr, axis=0)
        return b_arr_avg

    def get_mle(self, x):
        samples = np.random.exponential(self.params_tau_seed, size=len(x))
        return np.mean(samples)

    def fit_exponential(self, x, y):
        p = [
            self.params_mean_ion / 10,
            self.params_tau_seed / 5,
            self.params_mean_bkgnd,
        ]
        popt, pcov = curve_fit(Simulator.exponential_function, x, y, p0=p)
        return popt, pcov

    def plot_time_and_fit(self, x, y, popt, title="", display_fit=True):
        fig = plt.figure()
        ax = fig.gca()
        ax.step(x, y, label=title, where="post")
        yfit = Simulator.exponential_function(x, *popt)
        if display_fit:
            ax.plot(x, yfit, label="fit", alpha=0.6, color="#DC143C")  # color Crimson
        
        outfilename = f"{self.settings_output_path}{title}_ts{self.params_tau_seed:.2e}_t{popt[1]:.2e}"

        if self.settings_plot_titles:
            title = r'$\tau_{seed} =$' + f"{self.params_tau_seed:0.2e}" + ' [s], ' + r'$\tau =$' + f'{popt[1]:0.2e}'
        else:
            title=None
        
        ax.set(
            xlabel="Time [s]",
            ylabel="Amplitude [a.u.]",
            title=title,
        )
        ax.grid()
        ax.legend(loc="upper right", shadow=False)
        plt.tight_layout()
        plt.savefig(outfilename + ".png")
        plt.close()
        if self.settings_save_npz:
            np.savez(outfilename + ".npz", x=x, y=y, yfit=yfit)

    def fit_gauss_to_result_array(self, xvals, yvals, mu_est, sigma_est, title=""):
        hist = np.histogram(
            yvals,
            bins=1000,
            range=(-2 * self.params_tau_seed, 2 * self.params_tau_seed),
        )[0]

        print(hist)
        A_initial = max(hist)
        mu_initial = xvals[np.argmax(hist)]
        sigma_initial = np.std(hist)

        initial_guess = [A_initial, mu_initial, sigma_initial]
        print(initial_guess)

        # p = [1, mu_est, sigma_est]

        try:
            popt, pcov = curve_fit(
                Simulator.gaussian_function, xvals, hist, p0=initial_guess
            )
            yfit = Simulator.gaussian_function(xvals, *popt)
            fig = plt.figure()
            ax = fig.gca()
            ax.step(xvals, hist, label=title, where="post")

            ax.plot(xvals, yfit, label="fit", alpha=0.6, color="#DC143C")
            
            if self.settings_plot_titles:
                title=r'$\tau - \tau_{seed}$ [s]'
            else:
                title = None
            
            ax.set(
                xlabel=f"amp = {popt[0]:.2e}, mean = {popt[1]:.2e}, sigma = {abs(popt[2]):.2e}",
                ylabel="",
                title=title
            )

            ax.legend(loc="upper right", shadow=False)
            ax.grid()
            outfilename = f"{self.settings_output_path}{title}"
            plt.tight_layout()
            plt.savefig(outfilename + ".png")
            plt.close()
            if self.settings_save_npz:
                np.savez(outfilename + ".npz", xvals=xvals, y=hist, yfit=yfit)

        except Exception as e:
            logger.warning(f"{e}. Ignoring...")
            self.simulation_error_flag = True

    def fit_and_plot_gaussian(self, xvals, yvals, bins=100, title=""):
        counts, bin_edges = np.histogram(yvals, bins=bins)  # Binning y values
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute bin centers

        # Initial parameter estimates
        A_initial = max(counts)  # Peak amplitude
        mu_initial = bin_centers[np.argmax(counts)]  # Peak position
        sigma_initial = (
            max(bin_centers) - min(bin_centers)
        ) / 4  # Rough width estimate
        initial_guess = [A_initial, mu_initial, sigma_initial]

        try:
            # Fit Gaussian curve
            popt, _ = curve_fit(
                Simulator.gaussian_function, bin_centers, counts, p0=initial_guess
            )

            # Generate fitted curve
            x_fit = np.linspace(min(bin_centers), max(bin_centers), 100)
            y_fit = Simulator.gaussian_function(x_fit, *popt)

            # Plot results
            plt.figure(figsize=(8, 5))
            plt.hist(yvals, bins=bins, alpha=0.6, label=title)
            plt.plot(x_fit, y_fit, color="red", label="Fit", linewidth=2)
            # plt.xlabel("X values")
            plt.xlabel(
                f"amp = {popt[0]:.2e}, mean = {popt[1]:.2e}, sigma = {abs(popt[2]):.2e}"
            )
            if self.settings_plot_titles:
                plt.title("Gaussian Fit")
                
            plt.legend()

            # Save plot
            outfilename = f"{self.settings_output_path}{title}"
            plt.tight_layout()
            plt.savefig(outfilename + ".png")
            plt.close()

            if self.settings_save_npz:
                np.savez(outfilename + ".npz", x_fit=x_fit, y_fit=y_fit, yvals=yvals)

        except Exception as e:
            logger.warning(f"{e}. Ignoring...")
            self.simulation_error_flag = True

        return popt  # Returns (A, mu, sigma)

    def start(self):
        logger.info(f'Simulation tasks: {self.settings_tasks}')
        tau_events_arr = np.zeros(self.params_n_trials)
        sigma_events_arr = np.zeros(self.params_n_trials)

        tau_distro_arr = np.zeros(self.params_n_trials)
        sigma_distro_arr = np.zeros(self.params_n_trials)

        tau_mle_arr = np.zeros(self.params_n_trials)

        x = np.arange(
            0, self.params_timestep * self.params_n_sim_steps, self.params_timestep
        )  # start, stop, step in seconds

        xvals = np.arange(
            -self.params_tau_seed, self.params_tau_seed, 2 * self.params_tau_seed / 1000
        )

        if 'addup' in self.settings_tasks:
            logger.info('Performing simulation task: addup')
            for trial_number in tqdm(range(self.params_n_trials)):
                try:
                    # y1 = self.create_from_events_boxcar(x)
                    y1 = self.create_from_events_with_fluctuations(
                        x, trial_number, add_empty_shots=self.params_empty_shots
                    )
                    popt_events, pcov_events = self.fit_exponential(x, y1)
                    tau_events_arr[trial_number] = popt_events[1]

                    sigma_events_arr[trial_number] = np.sqrt(
                        pcov_events[1, 1]
                    )  # get error of taus from events

                    self.plot_time_and_fit(
                        x,
                        y1,
                        popt_events,
                        title=f"trial-{trial_number:04}_from_events",
                        display_fit=True,
                    )

                except (FloatingPointError, OptimizeWarning) as e:
                    logger.warning(e)
                    continue
                    
            logger.info("Creating tau_tru distribution from add up spectra.")
            popt = self.fit_and_plot_gaussian(
                xvals, self.params_tau_seed - tau_events_arr, title="tau_tru_from_addup"
            )
            sigma_tru_from_addup = popt[2]
            logger.info(f"Sigma_tru from add up is = {sigma_tru_from_addup}")

            logger.info("Creating sigma distribution from add up spectra.")
            popt = self.fit_and_plot_gaussian(
                xvals, sigma_events_arr, title="sigma_from_addup"
            )
            mean_of_sigmas_from_addup = popt[1]
            logger.info(f"Mean of sigmas from addup is = {mean_of_sigmas_from_addup}")

            rho_from_addup = mean_of_sigmas_from_addup / sigma_tru_from_addup
            logger.info(f"rho from addup = {rho_from_addup}")

            logger.info("Creating rho distribion from add up.")
            self.fit_and_plot_gaussian(
                xvals, sigma_events_arr / sigma_tru_from_addup, title="rho_from_addup"
            )

        if 'distro' in self.settings_tasks:
            logger.info('Performing simulation task: distro')
            for trial_number in tqdm(range(self.params_n_trials)):
                try:
                    y2 = self.create_from_distribution(x)
                    popt_distro, pcov_distro = self.fit_exponential(x, y2)
                    tau_distro_arr[trial_number] = popt_distro[1]
                    sigma_distro_arr[trial_number] = np.sqrt(pcov_distro[1,1]) # get error of taus from distro

                    self.plot_time_and_fit(
                        x, y2, popt_distro, title=f'trial-{trial_number:04}_from_distro', display_fit=True)

                except (FloatingPointError, OptimizeWarning) as e:
                    logger.warning(e)
                    continue

            logger.info('Creating tau_tru from binned exponential distribution.')
            popt = self.fit_and_plot_gaussian(xvals, self.params_tau_seed - tau_distro_arr, title = 'tau_tru_from_distro')
            sigma_tru_from_distro = popt[2]
            logger.info(f'Sigma_tru from distro is = {sigma_tru_from_distro}')

            logger.info('Creating sigma distribution from binned exponential distribution.')
            popt = self.fit_and_plot_gaussian(xvals, sigma_distro_arr, title='sigma_from_distro')
            mean_of_sigmas_from_distro = popt[1]
            logger.info(f'Mean of sigmas from distro = {mean_of_sigmas_from_distro}')

            rho_from_distro = mean_of_sigmas_from_distro / sigma_tru_from_distro
            logger.info(f'rho from distro = {rho_from_distro}')

            logger.info('Creating rho distribion from distro.')
            self.fit_and_plot_gaussian(xvals, sigma_distro_arr / sigma_tru_from_distro, title = 'rho_from_distro')

        if 'mle' in self.settings_tasks:
            logger.info('Performing simulation task: mle')
            for trial_number in tqdm(range(self.params_n_trials)):
                try:
                    tau_mle_arr[trial_number] = self.get_mle(x)
                    
                except (FloatingPointError, OptimizeWarning) as e:
                    logger.warning(e)
                    continue

            logger.info("Creating tau_tru from unbinned MLE distribion.")
            self.fit_and_plot_gaussian(
            xvals, self.params_tau_seed - tau_mle_arr, title="tau_tru_from_mle"
            )
    

        # Scatter plot
        logger.info("Creating the scatter plot.")
        fig, axs = plt.subplots()
        if 'addup' in self.settings_tasks:
            axs.step(
                np.arange(len(tau_events_arr)),
                self.params_tau_seed - tau_events_arr,
                where="post",
                label="tau_from_trials",
            )
        if 'distro' in self.settings_tasks:
            axs.step(np.arange(len(tau_distro_arr)), self.params_tau_seed -
                    tau_distro_arr, where='post', label='tau_from_distro')
        if 'mle' in self.settings_tasks:
            axs.step(
                np.arange(len(tau_mle_arr)),
                self.params_tau_seed - tau_mle_arr,
                where="post",
                label="tau_from_mle",
            )
        axs.set(
            xlabel="No. of trials",
            ylabel=r"$\tau_{seed} - \tau_{sim}$",
            title=r"$\tau_{seed} =$" + f"{self.params_tau_seed:0.2e}" + " [s]",
        )
        axs.legend(loc="upper right", shadow=False)
        plt.tight_layout()
        plt.savefig(f"{self.settings_output_path}scatter.png")
        plt.close()

        if self.simulation_error_flag:
            logger.warning(
                "There were some errors during the simulation. You might need to start again."
            )

# -----------------

class Params(BaseModel):
    tau_seed: float
    timestep: float
    n_sim_steps: int
    n_sim: int
    n_trials: int
    n_decay_steps: int = Field(..., ge=0, le=1)
    mean_ion: float
    stdv_ion: float
    mean_bkgnd: float
    stdv_bkgnd: float
    empty_shots: bool
    empty_shots_probability: float = Field(..., ge=0.0, le=1.0)


class Settings(BaseModel):
    output_path: str
    plot_every_event: bool
    save_npz: bool
    tasks: list
    plot_titles: bool


class Config(BaseModel):
    params: Params
    settings: Settings


def load_and_validate_toml(file_path: str) -> Config:
    data = toml.load(file_path)
    return Config(**data)

# -----------------

def main():
    logger.remove(0)
    logger.add(sys.stdout, level="INFO")

    parser = argparse.ArgumentParser(prog="decaysim")
    parser.add_argument(
        "param",
        nargs=1,
        type=str,
        default=None,
        help="Path and name of the simulation parameter file.",
    )

    # read command line args
    args = parser.parse_args()

    # load and check parameter file
    try:
        config = load_and_validate_toml(args.param[0])

    
    except ValueError as e:
        logger.error("Parameter file: " + str(e))
        exit()

    simulator = Simulator(config)
    try:
        simulator.start()
    except (EOFError, KeyboardInterrupt):
        logger.success("\nUser input cancelled. Aborting...")

    logger.info("Simulation end.")


# -----------------

if __name__ == "__main__":
    main()

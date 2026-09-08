#!/usr/bin/env python
"""
A particle decay time simulator for the investigation of the spectral add up method 

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
        self.simulation_duration = self.params_timestep * self.params_n_sim_steps
        
        self.params_n_sim = config.params.n_sim
        self.params_n_trials = config.params.n_trials
        self.params_n_decay_steps = config.params.n_decay_steps
        self.params_mean_ion_number = config.params.mean_ion_number
        self.params_max_ions = config.params.max_ions
        self.params_mean_ion = config.params.mean_ion
        self.params_noise_model = config.params.noise_model
        self.params_stdv_ion = config.params.stdv_ion
        self.params_mean_bkgnd = config.params.mean_bkgnd
        self.params_stdv_bkgnd = config.params.stdv_bkgnd
        self.params_empty_shots = config.params.empty_shots
        self.params_empty_shots_probability = config.params.empty_shots_probability
        self.params_limit = config.params.limit

        # Explicit, seedable RNG (no hidden global np.random state). A
        # fixed integer seed makes every run bit-reproducible; None falls
        # back to nondeterministic OS entropy, same as before this fix.
        self.params_seed = config.params.seed
        self.rng = np.random.default_rng(self.params_seed)

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
        logger.info(f"RNG seed: {self.params_seed if self.params_seed is not None else 'None (nondeterministic)'}")

    @staticmethod
    def gaussian_function(x, *p):
        return p[0] * np.exp(-((x - p[1]) ** 2) / (2.0 * p[2] ** 2))

    # @staticmethod
    # def exponential_function(x, *p):
    #     # Exponential function + offset
    #     return p[0] * np.exp(-x / p[1]) + p[2]  # amplitude, tau and constant

    # @staticmethod
    # def expo_func(x, A, B):
    #     return A * np.exp(-B * x)

    @staticmethod
    def expo_func_2(x, A, B):
        return A * np.exp(-x / B)

    @staticmethod
    def expo_func_3(x, A, B, C):
        return A * np.exp(-x / B) + C

    def create_from_events_boxcar(self, x):
        b_arr = np.array([])

        for i in range(self.params_n_sim):
            # create a single random number based on an exponential distribution
            num = self.rng.exponential(self.params_tau_seed)

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
        data = self.rng.exponential(self.params_tau_seed, size=len(x))
        hist = np.histogram(
            data,
            bins=int(
                self.params_timestep * self.params_n_sim_steps / self.params_timestep
            ),

            range=(0, self.params_timestep * self.params_n_sim_steps),

        )[0]
        
        return hist / np.max(hist)
    
    def plot_multi_ion_decay(self, x, n_alive, I_t, b, trial_number, limit):

        # Identify plateau boundaries (where n_alive changes)
        changes = np.where(np.diff(n_alive) != 0)[0]
        boundaries = np.concatenate(([0], changes + 1, [len(n_alive)]))
        x_ms = x*1000
        timestep_ms = x_ms[1] - x_ms[0]   # your sampling step in ms
        last_s = boundaries[-2]
        last_e = boundaries[-1]
        last_plateau = b[last_s:last_e]

        mean_last = last_plateau.mean()
        std_last  = last_plateau.std(ddof=1)
        while x_ms[-1] < limit:
            x_ms = np.append(x_ms, x_ms[-1] + timestep_ms)
            new_point = np.random.normal(mean_last, std_last)
            b = np.append(b, new_point)
            n_alive = np.append(n_alive, n_alive[-1])  # keep last plateau
        fig, ax = plt.subplots(figsize=(10, 5))

        # Plot raw intensity
        ax.plot(x_ms, b, label='Intensity')
        mean_label_used = False
        # Plot plateau means + std bands
        for i in range(len(boundaries) - 1):
            s = boundaries[i]
            e = boundaries[i + 1]
            
            n = n_alive[s]  # ion count in this plateau
            plateau = b[s:e]
            if x_ms[e-1] < 0.001:   # skip first 1 ms
                continue

            mean = plateau.mean()
            std = plateau.std(ddof=1)

            # horizontal line for plotting the mean
            ax.hlines(mean, x_ms[s] - (x_ms[1]-x_ms[0]), x_ms[e-1], alpha=0.6, colors='red', linestyle="--", label="Mean Intensity" if not mean_label_used else None)
            
            mean_label_used = True
            
            # annotate ion count
            interval_width = x_ms[e-1] - x_ms[s]
            
            if x_ms[e-1]>=limit:  # if the interval exceeds the x limit showed, adjust the label accordingly.
                last_x = limit
            else: last_x = x_ms[e-1]     #else:  Normal placement: centered in plateau
            label_x = (x_ms[s] + last_x) / 2
            #label_y = mean + 0.6*std
            ax.text(label_x, np.max(b),
                    f"{n} ions", ha='center', va='bottom',
                    fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
        
        while x_ms[-1] < limit:
            b = np.append(b, b[-1])  # or baseline value
            x_ms = np.append(x_ms, x_ms[-1] + (x_ms[1] - x_ms[0]))
        # Mark decay times
        decay_times = x_ms[changes]
        for t in decay_times:
            ax.axvline(t, color='grey', linestyle='--', linewidth=1)

        #ax.set_title(f"Multi-ion decay (trial {trial_number})")
        ax.set_xlim(0,limit)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Amplitude [a.u.]")
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax.set_xticks(plt.xticks()[0], plt.xticks()[0] / 1000.0)
        ax.yaxis.get_offset_text().set_fontsize(12)
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        outname = f"{self.settings_output_path}/multi_decay_trial_{trial_number:04}.png"
        plt.savefig(outname, dpi=150)
        plt.close()


    def create_from_events_with_fluctuations(
        self, x, trial_number, add_empty_shots=False
    ):
        """Simulate one add-up trial (n_sim injections, summed).

        noise_model options:
          "sqrt"       sigma_ion = stdv_ion * sqrt(n)
          "linear"     sigma_ion = stdv_bkgnd + n * stdv_ion
          "quadrature" sigma_ion = sqrt(stdv_bkgnd**2 + (n * stdv_ion)**2)
                       (RionSiS noise model, originally established in
                       Moritz Porstendoerfer's Bachelor thesis)
        """
        b_arr = np.array([])
        empty_shot_mask = self.rng.choice(
            [0, 1],
            size=self.params_n_sim,
            p=[
                self.params_empty_shots_probability,
                1 - self.params_empty_shots_probability,
            ],
        )

        t_max = 10 * self.params_tau_seed
        mask = x <= t_max

        for i in range(self.params_n_sim):

            # -----------------------------
            # 1. Draw number of ions
            # -----------------------------
            # Option A: Poisson (default)
            n_ions = self.rng.poisson(self.params_mean_ion_number)
            n_ions = min(n_ions, self.params_max_ions)


            # If empty shot → force zero ions
            if add_empty_shots and empty_shot_mask[i]:
                n_ions = 0

            # -----------------------------
            # 2. If no ions → pure background
            # -----------------------------
            if n_ions == 0:
                b = self.rng.normal(self.params_mean_bkgnd,
                                    self.params_stdv_bkgnd,
                                    len(x))
                #b = b[mask]
                b_arr = np.append(b_arr, b)
                continue

            decay_times = self.rng.exponential(self.params_tau_seed, size=n_ions)
            decay_times.sort()
            
            n_alive = n_ions - np.searchsorted(decay_times, x, side="right")
            
            # -----------------------------
            # 5. Convert n(t) to intensity using calibrated ladder
            # -----------------------------
            # I(n) = I0 + n*(I1 - I0)
            I0 = self.params_mean_bkgnd
            I1 = self.params_mean_ion
            I_t = I0 + n_alive * (I1 - I0)
            # other method:
            #I0 = self.rng.normal(self.params_mean_bkgnd, self.params_stdv_bkgnd)
            #I_t = n_alive * self.params_mean_ion + I0

            # ---------------------------------------------------------
            # Intermediate frame logic (multi-ion version)
            # ---------------------------------------------------------
            if self.params_n_decay_steps == 1:
                # For each ion, compute fractional alive time inside its decay bin
                for t_decay in decay_times:
                    frame = int(t_decay / self.params_timestep)
                    if frame < len(x):
                        # fractional time alive inside the bin
                        remainder = t_decay % self.params_timestep
                        fraction_alive = 1 - remainder / self.params_timestep

                        # subtract fractional amplitude from that bin
                        # but only if the ion is alive in that bin
                        if n_alive[frame] > 0:
                            I_t[frame] -= fraction_alive * (self.params_mean_ion - self.params_mean_bkgnd)

            
            # -----------------------------
            # 6. Noise model (selectable)
            # -----------------------------
            ladder_spacing = (I1 - I0)

            if self.params_noise_model == "sqrt":
                sigma_ion = self.params_stdv_ion * np.sqrt(n_alive)
                sigma_bkg = self.params_stdv_bkgnd

            elif self.params_noise_model == "linear":
                # sigma(n) = std_dev_bkg + n * std_dev_ion
                sigma_ion = self.params_stdv_bkgnd + n_alive * self.params_stdv_ion
                sigma_bkg = self.params_stdv_bkgnd

            elif self.params_noise_model == "quadrature":
                # sigma(n) = sqrt(std_dev_bkg^2 + (n * std_dev_ion)^2), i.e.
                # the RionSiS noise model (see docstring for provenance).
                sigma_ion = np.sqrt(self.params_stdv_bkgnd**2 + (n_alive * self.params_stdv_ion) ** 2)
                sigma_bkg = self.params_stdv_bkgnd

            else:
                # fallback: old model
                sigma_ion = self.params_stdv_ion * np.sqrt(n_alive)
                sigma_bkg = self.params_stdv_bkgnd

            # Zero-mean noise for ions and background, selected per bin by ion count
            noise_ion = self.rng.normal(0, sigma_ion)
            noise_bkg = self.rng.normal(0, sigma_bkg, len(x))
            noise = np.where(n_alive > 0, noise_ion, noise_bkg)

            # Add noise to ladder
            b = I_t + noise

            # Enforce positivity
            b = np.clip(b, I0 * 0.5, None)

            # Enforce ordering: I(n) > I(n-1)
            for n in range(1, np.max(n_alive) + 1):
                mask_n = (n_alive == n)
                mask_prev = (n_alive == n - 1)
                if np.any(mask_prev):
                    mean_prev = np.mean(b[mask_prev])
                    mean_curr = np.mean(b[mask_n])

                    if mean_curr <= mean_prev:
                    # shift entire plateau upward
                        shift = (mean_prev - mean_curr) + 0.1 * ladder_spacing
                        b[mask_n] += shift
            b_arr = np.append(b_arr, b)

            if self.settings_plot_every_event:
                x = np.arange(len(b)) * self.params_timestep
                y = b
                plt.step(x, y, where="post")  # color teal

                if self.settings_plot_titles:
                    plt.title(r'$\tau_{seed} =$'+str(self.params_tau_seed) + ' [s]')

                plt.xlabel("Time [s]")
                plt.ylabel("Amplitude [a.u.]")
                plt.xticks(plt.xticks()[0], plt.xticks()[0] / 1000.0)
                #plt.grid()
                outfilename = f"{self.settings_output_path}trial{trial_number:04}_decay{i:04}_ts{self.params_tau_seed:.2e}"
                plt.tight_layout()
                plt.savefig(outfilename + ".png")
                plt.close()
                if self.settings_save_npz:
                    np.savez(outfilename + ".npz", x=x, y=y)
        

            if trial_number < 50 and i==0:
                self.plot_multi_ion_decay(x[mask], n_alive[mask], I_t[mask], b[mask], trial_number, self.params_limit)

        b_arr = np.reshape(b_arr, (self.params_n_sim, len(x)))
        b_arr_avg = np.average(b_arr, axis=0)
        b_arr_sum = np.sum(b_arr, axis=0)
        return b_arr_sum

    def get_mle(self, x):
        # MLE of the exponential mean under Type-I right censoring at
        # obs_window: samples >= obs_window are censored (only known to
        # have survived at least obs_window), which the sample mean of
        # observed-only points ignores, biasing the estimate low.
        samples = self.rng.exponential(self.params_tau_seed, size=len(x))
        obs_window = self.params_timestep * self.params_n_sim_steps
        cutoff_samples = samples[samples < obs_window]
        n_censored = len(samples) - len(cutoff_samples)
        mle_est = np.mean(cutoff_samples) + n_censored * obs_window / len(cutoff_samples)
        return mle_est

    def fit_exponential(self, x, y):
        p = [
            #self.params_mean_ion,
            y[0],
            self.params_tau_seed,
            self.params_mean_bkgnd,
        ]
        
        #p = [1e-6, 0.1] 
        #popt, pcov = curve_fit(Simulator.exponential_function, x, y, p0=p)
        popt, pcov = curve_fit(Simulator.expo_func_3, x, y, p0=p)

        # Debug output
        # print('p_in: ', p)
        # print('p_out: ', popt)
        # print("Tau (τ) =", popt[1], "±", np.sqrt(pcov[1, 1]) / (popt[1]**2))
        # if popt[1] > 100:
        #     logger.error('Too large tau!')
        return popt, pcov

    def plot_time_and_fit(self, x, y, popt, limit, id_string="", display_fit=True):
        fig = plt.figure()
        ax = fig.gca()
        ax.step(x, y, label=id_string, where="post")
        yfit = Simulator.expo_func_3(x, *popt)
        if display_fit:
            ax.plot(x, yfit, label="fit",
                    #fr'$\tau =$ {popt[1]:0.2e}',
                    alpha=0.6, color="red", linestyle="--")  # color Crimson
        
        outfilename = f"{self.settings_output_path}{id_string}_ts{self.params_tau_seed:.2e}_t{popt[1]:.2e}_{self.simulation_duration:.2f}s"

        if self.settings_plot_titles:
            title = r'$\tau_{seed} =$' + f"{self.params_tau_seed:0.2e}" + ' [s], ' + r'$\tau =$' + f'{popt[1]:0.2e}'
        else:
            title=None
        
        ax.set(
            xlabel="Time [s]",
            ylabel="Amplitude [a.u.]",
            title=title,
        )
        # ax.grid()
        ax.set_xlim(-0.02,limit/1000)
        ax.grid(True)
        legend = ax.legend(fontsize=14)
        for text in legend.get_texts():
            text.set_verticalalignment('center')  # Options: 'top', 'bottom', 'center', 'baseline'
        
        plt.tight_layout()
        plt.savefig(outfilename + ".png")
        plt.close()
        if self.settings_save_npz:
            np.savez(outfilename + ".npz", x=x, y=y, yfit=yfit)

    def fit_and_plot_gaussian(self, xvals, yvals, bins=100, id_string=""):
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
            fig = plt.figure(figsize=(8, 5))
            ax = fig.gca()
            ax.hist(yvals, bins=bins, alpha=0.6, label=id_string)
            fit_line_label = fr'$\mu =$ {popt[1]:0.2e}' + "\n" + fr'$\sigma =$ {abs(popt[2]):.2e}'
            ax.plot(x_fit, y_fit, color="red", label=fit_line_label, linewidth=2)

            if self.settings_plot_titles:
                plt.title("Gaussian Fit")
                plt.xlabel(
                    f"amp = {popt[0]:.2e}, mean = {popt[1]:.2e}, sigma = {abs(popt[2]):.2e}"
                )
                
            legend = ax.legend(fontsize=14)
            for text in legend.get_texts():
                text.set_verticalalignment('center')  # Options: 'top', 'bottom', 'center', 'baseline'

            # Save plot
            outfilename = f"{self.settings_output_path}{id_string}_{self.simulation_duration:.2f}s"
            plt.tight_layout()
            plt.savefig(outfilename + ".png", dpi=300)
            plt.savefig(outfilename + ".svg")
            plt.close()

            if self.settings_save_npz:
                np.savez(outfilename + ".npz", x_fit=x_fit, y_fit=y_fit, yvals=yvals)

        except Exception as e:
            logger.error(f"{e}. Can't continue any more. Please rerun the simulation.")
            sys.exit(1)

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
        
        t_max = 10 * self.params_tau_seed
        mask = x <= t_max
        x_trunc = x[mask]
        
        xvals = np.arange(
            -self.params_tau_seed, self.params_tau_seed, 2 * self.params_tau_seed / 1000
        )

        if 'addup' in self.settings_tasks:
            logger.info('Performing simulation task: addup')
            for trial_number in tqdm(range(self.params_n_trials)):
            # Diagnostic: inspect tau_events_arr and sigma_events_arr
                

                try:
                    #y1 = self.create_from_events_boxcar(x)
                    y1 = self.create_from_events_with_fluctuations(
                        x, trial_number, add_empty_shots=self.params_empty_shots
                    )
                    x_trunc = x[x <= 10 * self.params_tau_seed]
                    popt_events, pcov_events = self.fit_exponential(x, y1)
                    tau_events_arr[trial_number] = popt_events[1]
                    
                    sigma_events_arr[trial_number] = np.sqrt(
                        pcov_events[1, 1]
                    )  # get error of taus from events

                    self.plot_time_and_fit(
                        x,
                        y1,
                        popt_events,
                        limit=self.params_limit,
                        id_string=f"trial-{trial_number:04}_from_events",
                        display_fit=True
                    )

                except (FloatingPointError, OptimizeWarning) as e:
                    logger.warning(e)
                    continue
                    
            logger.info("Creating tau_tru distribution from add up spectra.")
            # Debug output
            # print(tau_events_arr)
            popt = self.fit_and_plot_gaussian(
                xvals, self.params_tau_seed - tau_events_arr, id_string="tau_tru_from_addup"
            )
            sigma_tru_from_addup = popt[2]
            logger.info(f"Sigma_tru from add up is = {sigma_tru_from_addup}")

            logger.info("Creating sigma distribution from add up spectra.")
            
            plt.figure()
            plt.hist(sigma_events_arr, bins=50, alpha=0.7)
            plt.xlabel("σ_sim_old")
            plt.ylabel("Counts")
            plt.title("Distribution of σ_sim_old")
            plt.tight_layout()
            plt.savefig(f"{self.settings_output_path}/sigma_sim_old_distribution.png", dpi=150)
            plt.close()

            popt = self.fit_and_plot_gaussian(
                xvals, sigma_events_arr, id_string="sigma_from_addup"
            )
            mean_of_sigmas_from_addup = popt[1]
            logger.info(f"Mean of sigmas from addup is = {mean_of_sigmas_from_addup}")

            # new definition
            rho_from_addup = sigma_tru_from_addup / mean_of_sigmas_from_addup
            logger.info(f"rho from addup = {rho_from_addup}")

            # new definition
            logger.info("Creating the P distribion from add up.")
            self.fit_and_plot_gaussian(
                xvals, (self.params_tau_seed - tau_events_arr) / sigma_events_arr , id_string="P_from_addup"
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
                        x, y2, popt_distro, limit, id_string=f'trial-{trial_number:04}_from_distro', display_fit=True)

                except (FloatingPointError, OptimizeWarning) as e:
                    logger.warning(e)
                    continue

            logger.info('Creating tau_tru from binned exponential distribution.')
            popt = self.fit_and_plot_gaussian(xvals, self.params_tau_seed - tau_distro_arr, id_string = 'tau_tru_from_distro')
            sigma_tru_from_distro = popt[2]
            logger.info(f'Sigma_tru from distro is = {sigma_tru_from_distro}')

            logger.info('Creating sigma distribution from binned exponential distribution.')
            popt = self.fit_and_plot_gaussian(xvals, sigma_distro_arr, id_string='sigma_from_distro')
            mean_of_sigmas_from_distro = popt[1]
            logger.info(f'Mean of sigmas from distro = {mean_of_sigmas_from_distro}')

            # new definition
            rho_from_distro =  sigma_tru_from_distro / mean_of_sigmas_from_distro
            logger.info(f'rho from distro = {rho_from_distro}')

            # new definition
            logger.info('Creating the P distribion from distro.')
            self.fit_and_plot_gaussian(xvals, (self.params_tau_seed - tau_distro_arr) / sigma_distro_arr, id_string = 'P_from_distro')

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
            xvals, self.params_tau_seed - tau_mle_arr, id_string="tau_tru_from_mle"
            )
    

        # Scatter plot
        logger.info("Creating the scatter plot.")
        fig, axs = plt.subplots()
        if 'addup' in self.settings_tasks:
            axs.step(
                np.arange(len(tau_events_arr)),
                np.abs(self.params_tau_seed - tau_events_arr),
                where="post",
                label="tau_from_trials",
            )
        if 'distro' in self.settings_tasks:
            axs.step(np.arange(len(tau_distro_arr)), np.abs(self.params_tau_seed -
                    tau_distro_arr), where='post', label='tau_from_distro')
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
        axs.legend(fontsize=14)
        plt.tight_layout()
        outfilename = f"{self.settings_output_path}scatter_{self.simulation_duration:.2f}s"
        plt.savefig(outfilename + ".png", dpi=300)
        plt.savefig(outfilename + ".svg")
        plt.close()

# -----------------

class Params(BaseModel):
    tau_seed: float
    timestep: float
    n_sim_steps: int
    n_sim: int
    n_trials: int
    n_decay_steps: int = Field(..., ge=0, le=1)
    mean_ion_number: int
    noise_model: str
    max_ions: int
    mean_ion: float
    stdv_ion: float
    mean_bkgnd: float
    stdv_bkgnd: float
    empty_shots: bool
    empty_shots_probability: float = Field(..., ge=0.0, le=1.0)
    limit: int
    seed: int | None = None

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

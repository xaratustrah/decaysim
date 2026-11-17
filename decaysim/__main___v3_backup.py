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
from matplotlib.transforms import ScaledTranslation
import warnings
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.signal import lfilter
import numpy as np
from pydantic import BaseModel, Field
from tqdm import tqdm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# ======================================================================================================
# Common figure format for the paper figures (plot_multi_ion_decay and plot_time_and_fit).
# FIG_WIDTH_IN must be the width the figure has in the LaTeX document: in Overleaf print \the\linewidth
# and divide the points by 72.27 (e.g. 345.0pt -> 4.77 in). The file is then included WITHOUT rescaling,
# so the font sizes below appear in the paper exactly as set.
# ======================================================================================================
FIG_WIDTH_IN = 6.0    # inches
FIG_HEIGHT_IN = 3.6   # inches
FONT_AXIS, FONT_TICKS, FONT_LEGEND, FONT_ANNOT = 11, 10, 10, 8   # pt
LINEWIDTH = 1.5
X_LEFT_S = -0.02      # s, left x-limit shared by both figures
DPI = 300
GAMMA = 1.3957        # Lorentz factor: lab -> centre-of-mass (rest) frame, tau_cm = tau_lab / GAMMA
plt.rcParams.update({
    "font.size": FONT_AXIS, "axes.labelsize": FONT_AXIS, "axes.titlesize": FONT_AXIS,
    "xtick.labelsize": FONT_TICKS, "ytick.labelsize": FONT_TICKS, "legend.fontsize": FONT_LEGEND,
    "lines.linewidth": LINEWIDTH, "savefig.dpi": DPI, "pdf.fonttype": 42,
})

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

        # Bins excluded from the exponential fit. Must be the same bins as in the fit of the real data
        # (the real fit drops the first and the third bin), otherwise sigma_est is not comparable.
        self.fit_drop = tuple(config.params.fit_drop)
        # Same estimator as the real fit: start time, weights and comparison windows
        self.fit_t_start = config.params.fit_t_start            # [s] first time used in the main fit
        self.fit_weights = config.params.fit_weights            # "none" | "model" | "exp"
        self.compare_t_starts = list(config.params.compare_t_starts)
        self.rms_windows = [(float(a), float(b) if b is not None else np.inf) for a, b in config.params.rms_windows]
        # optional time-correlated ion noise (AR(1)), sigma = n * stdv_ion_corr, correlation length corr_bins
        self.params_stdv_ion_corr = config.params.stdv_ion_corr
        self.params_corr_bins = config.params.corr_bins
        # bootstrap (see bootstrap_residuals / bootstrap_injections)
        self.boot_n = config.params.bootstrap_n
        self.boot_block = config.params.bootstrap_block
        self.boot_trials = config.params.bootstrap_trials
        self.exp_injections_file = config.params.exp_injections_file
        self._boot = {}                    # results: {(who, method, t0): [sigma_tau, ...]}
        # own generator: the bootstrap does not change the simulated data of the seeded run
        self.boot_rng = np.random.default_rng(None if config.params.seed is None else config.params.seed + 1000)
        self._weights = None
        self._real_popt = None

        # Explicit, seedable RNG (no hidden global np.random state). A
        # fixed integer seed makes every run bit-reproducible; None falls
        # back to nondeterministic OS entropy.
        self.params_seed = config.params.seed
        self.rng = np.random.default_rng(self.params_seed)

        # handle trailing slash and file path properly
        self.settings_output_path = os.path.join(config.settings.output_path, "")
        self.settings_output_tag = config.settings.output_tag   # prefix of every output file name (e.g. "long_")
        self.settings_figure_formats = list(config.settings.figure_formats)
        self.settings_plot_every_event = config.settings.plot_every_event
        self.settings_save_npz = config.settings.save_npz
        self.settings_tasks = config.settings.tasks
        self.settings_plot_titles = config.settings.plot_titles
        self.settings_n_trial_plots = config.settings.n_trial_plots

        self.params_ion_distribution = config.params.ion_distribution
        self.params_ion_distribution_file = config.params.ion_distribution_file
        self.empirical_ion_counts = None
        self._ion_idx = 0
        self._diagnostic_level_by_n = {}
        self.exp_decay_file = config.params.exp_decay_file
        self._res_stats = {}

        # Experimental summed curve. Preferred: .npz from the half-life notebook (y, sigma, first_bin);
        # old .npy (only y, bins after fit_drop) still works. Stored on the full simulation grid (NaN = no data).
        self.exp_decay_values, self.exp_y_full, self.exp_s_full = None, None, None
        if self.exp_decay_file:
            n = self.params_n_sim_steps
            if self.exp_decay_file.endswith(".npz"):
                d = np.load(self.exp_decay_file)
                y_e = np.asarray(d["y"], float)
                s_e = np.asarray(d["sigma"], float) if "sigma" in d.files else None
                idx = (np.asarray(d["bins"], int) if "bins" in d.files
                       else int(d["first_bin"]) + np.arange(len(y_e)))
                if "dt" in d.files and abs(float(d["dt"]) - self.params_timestep) > 1e-3 * self.params_timestep:
                    logger.warning(f"timestep {self.params_timestep} differs from the experimental dt {float(d['dt'])}")
                if "n_inj" in d.files and int(d["n_inj"]) != self.params_n_sim:
                    logger.warning(f"n_sim {self.params_n_sim} differs from the experimental n_inj {int(d['n_inj'])}")
            else:
                y_e, s_e = np.asarray(np.load(self.exp_decay_file), float), None
                idx = np.flatnonzero(self.fit_mask(n))[:len(y_e)]       # saved after the mask
            ok = idx < n
            self.exp_y_full = np.full(n, np.nan); self.exp_y_full[idx[ok]] = y_e[ok]
            if s_e is not None:
                self.exp_s_full = np.full(n, np.nan); self.exp_s_full[idx[ok]] = s_e[ok]
            self.exp_decay_values = y_e
        if self.fit_weights == "exp" and self.exp_s_full is None:
            logger.error("fit_weights = 'exp' needs an exp_decay_file (.npz) with 'sigma'.")
            exit()

        if self.params_ion_distribution == "empirical":
            if not self.params_ion_distribution_file:
                logger.error("ion_distribution='empirical' requires ion_distribution_file to be set.")
                exit()
            if not os.path.isfile(self.params_ion_distribution_file):
                logger.error(f"ion_distribution_file not found: {self.params_ion_distribution_file}")
                exit()
            self.empirical_ion_counts = np.load(self.params_ion_distribution_file)
            logger.info(
                f"Loaded empirical ion-count distribution: {len(self.empirical_ion_counts)} samples, "
                f"mean={self.empirical_ion_counts.mean():.3f}, max={int(self.empirical_ion_counts.max())}"
            )
        elif self.params_ion_distribution != "poisson":
            logger.error(f"Unknown ion_distribution: {self.params_ion_distribution!r} (use 'poisson' or 'empirical')")
            exit()

        if not os.path.isdir(self.settings_output_path):
            logger.error("Output file path does not exist.")
            exit()

        logger.info("Simulation start. Enter command or ctrl-C to abort.")
        logger.info(f"Output file path: {self.settings_output_path}")
        logger.info(f"RNG seed: {self.params_seed if self.params_seed is not None else 'None (nondeterministic)'}")

    # ------------------------------------------------------------------
    # static helpers
    # ------------------------------------------------------------------
    @staticmethod
    def gaussian_function(x, *p):
        return p[0] * np.exp(-((x - p[1]) ** 2) / (2.0 * p[2] ** 2))

    @staticmethod
    def expo_func_2(x, A, B):
        return A * np.exp(-x / B)

    @staticmethod
    def expo_func_3(x, A, B, C):
        return A * np.exp(-x / B) + C

    # ------------------------------------------------------------------
    # shared helpers: both paper figures are built and saved through these, so they cannot drift apart
    # ------------------------------------------------------------------
    def _outpath(self, stem):
        """Full path of an output file; the optional output_tag keeps files of different runs apart
        (important for \\includesvg, which caches converted files by base name only)."""
        return f"{self.settings_output_path}{self.settings_output_tag}{stem}"

    def _new_figure(self):
        return plt.subplots(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), constrained_layout=True)

    def _style_axes(self, ax, limit, legend_loc="upper right"):
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Amplitude [a.u.]")
        ax.set_xlim(X_LEFT_S, limit / 1000.0)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(FONT_TICKS)
        ax.grid(True)
        legend = ax.legend(fontsize=FONT_LEGEND, loc=legend_loc)
        for text in legend.get_texts():
            text.set_verticalalignment("center")

    def _save_figure(self, fig, outbase):
        for ext in self.settings_figure_formats:
            fig.savefig(f"{outbase}.{ext}", dpi=DPI)   # no bbox_inches="tight": exact figure size is kept
        plt.close(fig)

    def fit_mask(self, n):
        """True for the bins that enter the exponential fit (the bins in self.fit_drop are excluded)."""
        keep = np.ones(n, dtype=bool)
        keep[list(self.fit_drop)] = False
        return keep

    # ------------------------------------------------------------------
    # simulation of the decays
    # ------------------------------------------------------------------
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

    def draw_ion_count(self, rng: np.random.Generator) -> int:
        """Draw a single injection's initial ion count: Poisson, or next-in-order from the empirical array.
        max_ions is optional; if set, draws above it are CLIPPED to it (lowering the realized mean)."""
        if self.params_ion_distribution == "empirical":
            n_ions = int(self.empirical_ion_counts[self._ion_idx % len(self.empirical_ion_counts)])
            self._ion_idx += 1
        else:
            n_ions = int(rng.poisson(self.params_mean_ion_number))
        if self.params_max_ions is not None:
            n_ions = min(n_ions, self.params_max_ions)
        return n_ions

    # ------------------------------------------------------------------
    # figure 1: one simulated injection with the model staircase
    # ------------------------------------------------------------------
    def plot_multi_ion_decay(self, x, n_alive, I_t, b, sim_number, limit):
        x = np.asarray(x, dtype=float)
        b = np.asarray(b, dtype=float)
        n_alive = np.asarray(n_alive)
        limit_s = limit / 1000.0
        timestep = x[1] - x[0]
        changes = np.where(np.diff(n_alive) != 0)[0]
        boundaries = np.concatenate(([0], changes + 1, [len(n_alive)]))
        decay_x = x[changes]                      # before padding, so the indices stay valid

        # pad to the plotted range with the model noise of the last level. A LOCAL generator is used on purpose:
        # drawing from self.rng here would shift the seeded simulation stream (and change the results).
        pad_rng = np.random.default_rng(sim_number)
        n_last = n_alive[-1]
        pad_level = self.params_mean_bkgnd + n_last * (self.params_mean_ion - self.params_mean_bkgnd)
        pad_sigma = np.sqrt(self.params_stdv_bkgnd**2 + (n_last * self.params_stdv_ion) ** 2)
        while x[-1] < limit_s:
            x = np.append(x, x[-1] + timestep)
            b = np.append(b, pad_rng.normal(pad_level, pad_sigma))
            n_alive = np.append(n_alive, n_alive[-1])

        fig, ax = self._new_figure()
        ax.plot(x, b, color="C0", label="Simulated signal")
        level_label_used = False
        for i in range(len(boundaries) - 1):
            s, e = boundaries[i], boundaries[i + 1]
            if x[s] >= limit_s:
                continue
            n = n_alive[s]
            # model level I0 + n*(I1-I0): monotonic by construction (plateau MEANS of the noisy bins are not)
            level = self.params_mean_bkgnd + n * (self.params_mean_ion - self.params_mean_bkgnd)
            x_from = x[s] - timestep if s > 0 else x[s]
            x_to = x[-1] if i == len(boundaries) - 2 else max(x[e - 1], x_from + timestep)   # at least one bin wide
            ax.hlines(level, x_from, x_to, alpha=0.6, colors="red", linestyle="--",
                      label="Staircase level" if not level_label_used else None)
            level_label_used = True
            # label sits just above its own level line: the levels are distinct, so labels do not collide
            ax.text((x[s] + min(x_to, limit_s)) / 2, level, f"{n} ion{'s' if n != 1 else ''}",
                    transform=ax.transData + ScaledTranslation(0, 3 / 72, fig.dpi_scale_trans),
                    ha="center", va="bottom", fontsize=FONT_ANNOT, bbox=dict(facecolor="white", alpha=0.8, pad=1.0))
        for t in decay_x:
            ax.axvline(t, color="grey", linestyle="--", linewidth=1)

        self._style_axes(ax, limit, legend_loc="upper right")
        self._save_figure(fig, self._outpath(f"multi_decay_trial_{sim_number:04}"))

    def create_from_events_with_fluctuations(
        self, x, trial_number, add_empty_shots=False
    ):
        """Simulate one add-up trial (n_sim injections, summed).

        Returns (summed signal, sigma used for the fit, all decay times of this trial). The decay times are the
        real latent events the summed curve is built from; they are what the MLE benchmark must use.

        noise_model options:
          "sqrt"       sigma_ion = stdv_ion * sqrt(n)
          "linear"     sigma_ion = stdv_bkgnd + n * stdv_ion
          "quadrature" sigma_ion = sqrt(stdv_bkgnd**2 + (n * stdv_ion)**2)
                       (RionSiS noise model, originally established in
                       Moritz Porstendoerfer's Bachelor thesis)
        """
        if self.params_ion_distribution == "empirical":
            self._ion_idx = 0  # replay the real sequence in order, fresh for every trial

        b_arr = np.empty((self.params_n_sim, len(x)))      # preallocated (np.append in the loop was O(N^2))
        var_sum = np.zeros(len(x))                          # model variance of the SUM, bin by bin
        all_decay_times = []
        # mask == 0 with probability empty_shots_probability  ->  those injections are forced empty
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
            n_ions = self.draw_ion_count(self.rng)

            # If empty shot -> force zero ions
            if add_empty_shots and not empty_shot_mask[i]:
                n_ions = 0

            # -----------------------------
            # 2. If no ions -> pure background
            # -----------------------------
            if n_ions == 0:
                b = self.rng.normal(self.params_mean_bkgnd,
                                    self.params_stdv_bkgnd,
                                    len(x))
                b_arr[i] = b
                var_sum += self.params_stdv_bkgnd**2

                if trial_number == 1 and i < 50:
                    n_alive_empty = np.zeros(len(x), dtype=int)
                    I_t_empty = np.full(len(x), self.params_mean_bkgnd)
                    self.plot_multi_ion_decay(x[mask], n_alive_empty[mask], I_t_empty[mask], b[mask], i, self.params_limit)
                continue

            decay_times = self.rng.exponential(self.params_tau_seed, size=n_ions)
            decay_times.sort()
            all_decay_times.append(decay_times)

            n_alive = n_ions - np.searchsorted(decay_times, x, side="right")

            # -----------------------------
            # 3. Convert n(t) to intensity using calibrated ladder: I(n) = I0 + n*(I1 - I0)
            # -----------------------------
            I0 = self.params_mean_bkgnd
            I1 = self.params_mean_ion
            I_t = I0 + n_alive * (I1 - I0)

            # ---------------------------------------------------------
            # Intermediate frame logic (multi-ion version): the bin in which an ion decays only holds
            # the fraction of the bin during which it was alive (bin-averaged signal)
            # ---------------------------------------------------------
            if self.params_n_decay_steps == 1:
                for t_decay in decay_times:
                    frame = int(t_decay / self.params_timestep)
                    if frame < len(x):
                        remainder = t_decay % self.params_timestep
                        fraction_alive = 1 - remainder / self.params_timestep
                        if n_alive[frame] > 0:
                            I_t[frame] -= fraction_alive * (self.params_mean_ion - self.params_mean_bkgnd)

            # -----------------------------
            # 4. Noise model (selectable)
            # -----------------------------
            if self.params_noise_model == "sqrt":
                sigma_ion = self.params_stdv_ion * np.sqrt(n_alive)
                sigma_bkg = self.params_stdv_bkgnd

            elif self.params_noise_model == "linear":
                # sigma(n) = std_dev_bkg + n * std_dev_ion
                sigma_ion = self.params_stdv_bkgnd + n_alive * self.params_stdv_ion
                sigma_bkg = self.params_stdv_bkgnd

            elif self.params_noise_model == "quadrature":
                # sigma(n) = sqrt(std_dev_bkg^2 + (n * std_dev_ion)^2), i.e. the RionSiS noise model
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
            sigma_bin = np.where(n_alive > 0, sigma_ion, sigma_bkg)

            # optional time-correlated ion noise: AR(1) with unit variance, amplitude n * stdv_ion_corr
            if self.params_stdv_ion_corr > 0:
                a = np.exp(-1.0 / self.params_corr_bins) if self.params_corr_bins > 0 else 0.0
                w = self.rng.normal(0, 1, len(x))
                c = np.sqrt(1 - a * a)
                if c > 0:
                    w[0] /= c                                    # e[0] ~ N(0,1): stationary from the first bin
                e = lfilter([c], [1, -a], w) if c > 0 else np.full(len(x), w[0])   # e[j] = a e[j-1] + c w[j]
                noise = noise + n_alive * self.params_stdv_ion_corr * e
                sigma_bin = np.sqrt(sigma_bin**2 + (n_alive * self.params_stdv_ion_corr) ** 2)
            var_sum += sigma_bin**2

            # Add noise to ladder (no clipping and no forced ordering: both bias the add-up statistics)
            b = I_t + noise

            changes = np.where(np.diff(n_alive) != 0)[0]
            boundaries = np.concatenate(([0], changes + 1, [len(n_alive)]))
            for j in range(len(boundaries) - 1):
                s, e = boundaries[j], boundaries[j + 1]
                n = int(n_alive[s])
                plateau_mean = b[s:e].mean()
                self._diagnostic_level_by_n.setdefault(n, []).append(plateau_mean)
            b_arr[i] = b

            if self.settings_plot_every_event:
                x_ev = np.arange(len(b)) * self.params_timestep
                fig_ev, ax_ev = self._new_figure()
                ax_ev.step(x_ev, b, where="post")
                ax_ev.set_xlabel("Time [s]")
                ax_ev.set_ylabel("Amplitude [a.u.]")
                ax_ev.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
                ax_ev.yaxis.get_offset_text().set_fontsize(FONT_TICKS)
                if self.settings_plot_titles:
                    ax_ev.set_title(r"$\tau_{seed} =$" + str(self.params_tau_seed) + " [s]")
                outbase = self._outpath(f"trial{trial_number:04}_decay{i:04}_ts{self.params_tau_seed:.2e}")
                self._save_figure(fig_ev, outbase)
                if self.settings_save_npz:
                    np.savez(outbase + ".npz", x=x_ev, y=b)

            if trial_number == 1 and i < 50:
                self.plot_multi_ion_decay(x[mask], n_alive[mask], I_t[mask], b[mask], i, self.params_limit)

        b_arr_sum = np.sum(b_arr, axis=0)
        trial_decay_times = (
            np.concatenate(all_decay_times) if all_decay_times else np.array([])
        )
        # sigma of the SUMMED curve from the noise model (was: sigma of the last injection only)
        self._last_b_arr = b_arr                            # single injections of this trial (for the injection bootstrap)
        return b_arr_sum, np.sqrt(var_sum), trial_decay_times

    def get_mle(self, decay_times):
        """MLE of the exponential mean under Type-I right censoring at obs_window, computed from the ACTUAL
        decay times already simulated for this trial's add-up curve (not an independent redraw) -- this is what
        makes it a genuine same-data benchmark against sigma_tru from the add-up method."""
        if len(decay_times) == 0:
            logger.warning("get_mle: no decays this trial (all injections empty) -- returning nan.")
            return np.nan
        obs_window = self.params_timestep * self.params_n_sim_steps
        cutoff_samples = decay_times[decay_times < obs_window]
        n_censored = len(decay_times) - len(cutoff_samples)
        if len(cutoff_samples) == 0:
            logger.warning("get_mle: zero observed (uncensored) decays this trial -- returning nan.")
            return np.nan
        mle_est = np.mean(cutoff_samples) + n_censored * obs_window / len(cutoff_samples)
        return mle_est

    def fit_exponential(self, x, y, sigma=None, t_start=None):
        """Fit of A*exp(-t/tau)+C, identical to fit_exp() of the half-life notebook:
        bins in fit_drop and t < t_start excluded, same p0 logic and bounds, weighted if sigma is given
        (absolute_sigma=False -> errors scaled by sqrt(chi2_r), as in the notebook 'with σ' fit)."""
        t_start = self.fit_t_start if t_start is None else t_start
        keep = self.fit_mask(len(x)) & (x >= t_start - 1e-9) & np.isfinite(y)
        x, y = x[keep], y[keep]
        s = None if sigma is None else np.asarray(sigma)[keep]
        p0 = [max(y[0] - y[-1], 1e-15), self.params_tau_seed, y[-1]]
        popt, pcov = curve_fit(Simulator.expo_func_3, x, y, p0=p0, sigma=s, absolute_sigma=False,
                               bounds=([0, 1e-4, -np.inf], [np.inf, 50, np.inf]), maxfev=20000)
        logger.debug(f"params: A={popt[0]:.4e}, tau={popt[1]:.5f}, sigma_tau={np.sqrt(pcov[1, 1]):.5f}, C={popt[2]:.4e}")
        return popt, pcov

    def fit_stats(self, x, y, sigma=None, t_start=None):
        """tau, sigma_tau, chi2_r, lag-1 and residual rms in 3 time windows, for one fit (real or simulated)."""
        t_start = self.fit_t_start if t_start is None else t_start
        popt, pcov = self.fit_exponential(x, y, sigma, t_start)
        keep = self.fit_mask(len(x)) & (x >= t_start - 1e-9) & np.isfinite(y)
        t, r = x[keep], (y - Simulator.expo_func_3(x, *popt))[keep]
        s = np.ones_like(r) if sigma is None else np.asarray(sigma)[keep]
        rms = [np.sqrt(np.mean(r[(t >= lo) & (t < hi)] ** 2)) if np.any((t >= lo) & (t < hi)) else np.nan
               for lo, hi in self.rms_windows]
        chi2r = np.sum((r / s) ** 2) / (len(r) - 3) if sigma is not None else np.nan   # no meaning without σ
        return dict(tau=popt[1], dtau=np.sqrt(pcov[1, 1]), chi2r=chi2r,
                    lag1=np.corrcoef(r[:-1], r[1:])[0, 1], rms=rms)


    # ------------------------------------------------------------------
    # figure 2: summed-up signal, experimental data and fit
    # ------------------------------------------------------------------
    def plot_time_and_fit(self, x, y, popt, limit, id_string="", display_fit=True, show_excluded=False, pcov=None):
        keep = self.fit_mask(len(x))
        # main plot + residuals below (same width as the other figures, taller for the residual panel)
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN * 1.35), sharex=True,
                                     gridspec_kw=dict(height_ratios=[3, 1]), constrained_layout=True)
        # only the bins used in the fit (fit_drop skipped), at their real times: no shift towards t = 0
        ax.step(x[keep], y[keep], color="C0", label="Summed-up simulated signal", where="mid")
        if show_excluded:
            ax.plot(x[~keep], y[~keep], "o", mfc="none", color="C0", label="Excluded from fit")

        if self.exp_y_full is not None:
            ok = np.isfinite(self.exp_y_full)
            ax.plot(x[ok], self.exp_y_full[ok], color="C1", label="Experimental data")

        if display_fit:
            fit_label = "Simulated fit"
            if pcov is not None:                       # fitted lifetime (the value printed as "tau:") and rest-frame half-life
                tau_ms, dtau_ms = popt[1] * 1e3, np.sqrt(pcov[1, 1]) * 1e3
                fit_label += ("\n" + r"$\tau$" + f" = {tau_ms:.2f} ± {dtau_ms:.2f} ms"
                              + "\n" + r"$t_{1/2}^{\mathrm{rest}}$" + f" = {tau_ms*np.log(2)/GAMMA:.2f} ± {dtau_ms*np.log(2)/GAMMA:.2f} ms")
            tf = np.linspace(x[keep][0], x[-1], 400)
            ax.plot(tf, Simulator.expo_func_3(tf, *popt), label=fit_label, alpha=0.6, color="red", linestyle="--")

            # residuals (data − simulated fit) below the main plot
            res = y - Simulator.expo_func_3(x, *popt)
            keep = keep & (x >= self.fit_t_start - 1e-9)          # only the bins of the fit
            rx.step(x[keep], res[keep], color="C0", where="mid", label="Simulated − fit")
            if show_excluded:
                rx.plot(x[~keep], res[~keep], "o", mfc="none", color="C0")
            # residuals of the REAL data with the same fit (same bins, start and weights)
            if self.exp_y_full is not None:
                if self._real_popt is None:
                    try:
                        self._real_popt = self.fit_exponential(x, self.exp_y_full, self._weights)[0]
                    except Exception as e:
                        logger.warning(f"real fit for the residual plot failed: {e}")
                if self._real_popt is not None:
                    kr = keep & (x >= self.fit_t_start - 1e-9) & np.isfinite(self.exp_y_full)
                    res_r = self.exp_y_full - Simulator.expo_func_3(x, *self._real_popt)
                    rx.step(x[kr], res_r[kr], color="C1", where="mid", alpha=0.8, label="Experimental − fit")
            rx.axhline(0, color="red", linestyle="--", alpha=0.6)
            rx.legend(fontsize=FONT_ANNOT, loc="upper right", ncol=2)

        if self.settings_plot_titles:
            ax.set_title(r"$\tau_{seed} =$" + f"{self.params_tau_seed:0.2e}" + " [s], " + r"$\tau =$" + f"{popt[1]:0.2e}")

        self._style_axes(ax, limit, legend_loc="upper right")
        ymax = max(np.max(np.abs(y[keep])), np.nanmax(np.abs(self.exp_y_full)) if self.exp_y_full is not None else 0)
        order = 10 ** np.floor(np.log10(ymax))      # e.g. 1e-5 if the data peak at ~3e-5

        self._style_axes(ax, limit, legend_loc="upper right")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.yaxis.set_major_locator(MultipleLocator(order/2))   # ticks at 1, 2, 3 … ×order 
        ax.yaxis.get_major_formatter().set_useOffset(True)
        for text in ax.get_legend().get_texts():                    # multi-line fit label: keep it inside the legend box
            text.set_verticalalignment("baseline")
        ax.set_xlabel(""); ax.tick_params(labelbottom=False)        # time axis only on the residual panel
        rx.set_xlabel("Time [s]"); rx.set_ylabel("Residual [a.u.]")
        rx.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        rx.yaxis.get_offset_text().set_fontsize(FONT_TICKS)
        rx.grid(True)
        outbase = self._outpath(
            f"{id_string}_ts{self.params_tau_seed:.2e}_t{popt[1]:.2e}_{self.simulation_duration:.2f}s"
        )
        self._save_figure(fig, outbase)
        if self.settings_save_npz:
            np.savez(outbase + ".npz", x=x, y=y, yfit=Simulator.expo_func_3(x, *popt), fit_keep=keep)

    def fit_and_plot_gaussian(self, xvals, yvals, bins=100, id_string="", label=None,
                          xlabel="Value [a.u.]", ylabel="Counts"):
        label = label or id_string
        counts, bin_edges = np.histogram(yvals, bins=bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        initial_guess = [max(counts), bin_centers[np.argmax(counts)],
                        (max(bin_centers) - min(bin_centers)) / 4]

        try:
            popt, _ = curve_fit(Simulator.gaussian_function, bin_centers, counts, p0=initial_guess)
            A, mu, sigma = popt[0], popt[1], abs(popt[2])

            x_fit = np.linspace(min(bin_centers), max(bin_centers), 400)
            y_fit = Simulator.gaussian_function(x_fit, *popt)

            fig = plt.figure(figsize=(8, 5))
            ax = fig.gca()
            ax.hist(yvals, bins=bins, alpha=0.6, label=label)
            ax.plot(x_fit, y_fit, color="red", linewidth=2,
                    label="Gaussian fit\n" + fr"$\mu =$ {mu:0.2e}" + "\n" + fr"$\sigma =$ {sigma:.2e}")

            # dashed line at the centre of the distribution
            ax.axvline(mu, color="k", linestyle="--", linewidth=1.2)

            # sigma shown directly in the plot: arrow from mu to mu+sigma at the 1-sigma height
            y_sig = A * np.exp(-0.5)
            ax.annotate("", xy=(mu + sigma, y_sig), xytext=(mu, y_sig),
                        arrowprops=dict(arrowstyle="<->", color="k", lw=1.5))
            ax.text(mu + sigma / 2, y_sig, r"$\sigma$", ha="center", va="bottom", fontsize=16,
                    bbox=dict(facecolor="white", alpha=0.7, pad=1.0, edgecolor="none"))

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if self.settings_plot_titles:
                ax.set_title("Gaussian Fit")

            legend = ax.legend(fontsize=16, loc="upper right")                 # was 14
            #for text in legend.get_texts():
            #    text.set_verticalalignment('center')

            outfilename = self._outpath(f"{id_string}_{self.simulation_duration:.2f}s")
            plt.tight_layout()
            plt.savefig(outfilename + ".png", dpi=300)
            plt.savefig(outfilename + ".svg")
            plt.close()

            if self.settings_save_npz:
                np.savez(outfilename + ".npz", x_fit=x_fit, y_fit=y_fit, yvals=yvals)

        except Exception as e:
            # do not stop the program: fall back to the plain mean and std of the values
            yv = np.asarray(yvals, float); yv = yv[np.isfinite(yv)]
            logger.warning(f"Gaussian fit of {id_string} failed ({e}): using mean/std of the values instead.")
            popt = np.array([len(yv), np.mean(yv), np.std(yv)])

        return popt

    # ------------------------------------------------------------------
    # bootstrap: two ways to get the uncertainty of tau from ONE data set
    # ------------------------------------------------------------------
    def bootstrap_residuals(self, x, y, sigma, t0):
        """Residual (moving-block) bootstrap of the summed curve: y* = fit + resampled residuals, refit, std(tau*).
        Residuals are divided by sigma before resampling (noise changes with time) and multiplied back after.
        bootstrap_block = 1: plain bootstrap (independent bins); > 1: blocks of that many bins keep the
        bin-to-bin correlation of the residuals. Only needs the summed curve, but it only knows the noise
        that is LEFT in the residuals: fluctuations that look like a slightly different exponential
        (e.g. the random decay times) are absorbed by the fit and are NOT included."""
        popt, _ = self.fit_exponential(x, y, sigma, t0)
        keep = self.fit_mask(len(x)) & (x >= t0 - 1e-9) & np.isfinite(y)
        s = np.ones(keep.sum()) if sigma is None else np.asarray(sigma)[keep]
        f = Simulator.expo_func_3(x, *popt)
        z = (y[keep] - f[keep]) / s
        n, L = len(z), max(1, min(self.boot_block, len(z)))
        taus = []
        for _ in range(self.boot_n):
            starts = self.boot_rng.integers(0, n - L + 1, int(np.ceil(n / L)))
            zs = np.concatenate([z[k:k + L] for k in starts])[:n]
            ys = np.full(len(x), np.nan); ys[keep] = f[keep] + zs * s
            try:
                taus.append(self.fit_exponential(x, ys, sigma, t0)[0][1])
            except Exception:
                pass
        return np.std(taus, ddof=1) if len(taus) > 2 else np.nan

    def bootstrap_injections(self, x, Y, sigma, t0):
        """Injection bootstrap: draw the injections (rows of Y) with replacement, sum, refit, std(tau*).
        Includes EVERYTHING that changes from injection to injection (number of ions, decay times, noise),
        so it estimates the true spread directly. Needs the single-injection curves."""
        nI, taus = Y.shape[0], []
        for _ in range(self.boot_n):
            cnt = np.bincount(self.boot_rng.integers(0, nI, nI), minlength=nI).astype(float)
            ys = cnt @ Y
            try:
                taus.append(self.fit_exponential(x, ys, sigma, t0)[0][1])
            except Exception:
                pass
        return np.std(taus, ddof=1) if len(taus) > 2 else np.nan

    def load_exp_injections(self):
        """Single-injection curves of the REAL data: npz with Y [n_inj, n_bins] and bins (time-bin index of every
        column, same grid as the summed curve) or first_bin, and dt. Returns (x, Y on the full grid, NaN = no data)."""
        d = np.load(self.exp_injections_file)
        Y = np.asarray(d["Y"], float)
        bins = np.asarray(d["bins"], int) if "bins" in d.files else int(d["first_bin"]) + np.arange(Y.shape[1])
        dt = float(d["dt"]) if "dt" in d.files else self.params_timestep
        n = int(bins.max()) + 1
        Yf = np.full((Y.shape[0], n), np.nan); Yf[:, bins] = Y
        good = np.all(np.isfinite(Y), axis=1)
        logger.info(f"real single injections: {good.sum()} of {Y.shape[0]} used, {Y.shape[1]} bins, dt = {dt*1e3:.3f} ms")
        return np.arange(n) * dt, Yf[good]

    # ------------------------------------------------------------------
    # real vs simulated: same fit, same bins, same weights
    # ------------------------------------------------------------------
    def bootstrap_table(self, x, weights, tau_events_arr):
        """Uncertainty of tau_lab from every method, real and simulated, side by side."""
        B = lambda who, m, t0: np.nanmean(self._boot.get((who, m, t0), [np.nan]))
        L = (["=" * 96,
              f"UNCERTAINTY OF τ_lab [ms]: fit error vs bootstrap vs RionDS  ({self.boot_n} resamples, "
              f"sim bootstrap on the first {self.boot_trials} trials, block = {self.boot_block} bins)",
              "  residual bootstrap : resample the residuals of the summed curve (needs only the summed curve)",
              "  injection bootstrap: resample the injections and re-sum (needs the single-injection curves)",
              "  RionDS σ_true      : spread of τ over the simulated trials (true answer for the simulation)",
              "=" * 96])
        for t0 in self.compare_t_starts:
            S = self._res_stats.get(t0, [])
            tau_s = np.array([s_["tau"] for s_ in S]); dt_s = np.array([s_["dtau"] for s_ in S])
            sig_true = np.nanstd(tau_s)
            R = None
            if self.exp_y_full is not None:
                try:
                    R = self.fit_stats(x, self.exp_y_full, weights, t0)
                except Exception:
                    pass
            rows = [("fit error σ_fit", R["dtau"] if R else np.nan, np.nanmean(dt_s)),
                    (f"residual bootstrap (block {self.boot_block})", B("real", "res", t0), B("sim", "res", t0)),
                    ("injection bootstrap", B("real", "inj", t0), B("sim", "inj", t0)),
                    ("RionDS σ_true", np.nan, sig_true)]
            L.append(f"-- fit from {t0*1e3:.0f} ms --")
            L.append(f"  {'method':32}{'REAL':>12}{'SIM':>12}{'SIM / σ_true':>16}")
            for name, r_, s_ in rows:
                f_ = lambda v: f"{v*1e3:12.3f}" if np.isfinite(v) else f"{'–':>12}"
                L.append(f"  {name:32}{f_(r_)}{f_(s_)}{(s_/sig_true if np.isfinite(s_) else np.nan):16.2f}")
            if R:
                rho = sig_true / np.nanmean(dt_s)
                L.append(f"  -> real uncertainty from RionDS: σ_fit(real) × ρ = {R['dtau']*1e3:.3f} × {rho:.2f} = "
                         f"{R['dtau']*rho*1e3:.3f} ms")
            if self.exp_injections_file and self.exp_injections_file != self.exp_decay_file:
                L.append("  (real injection bootstrap uses the single-injection file: its own band/binning, see log)")
        L.append("=" * 96)
        return L

    def print_comparison(self, x, weights, tau_events_arr):
        w = {"none": "no σ", "model": "with σ (model)", "exp": "with σ (exp. σ)"}[self.fit_weights]
        lines = ["=" * 96,
                 f"REAL vs SIMULATED  (fit {w}, bins {list(self.fit_drop)} dropped; sim = mean ± std over trials)",
                 "  σ_fit = fit error of τ_lab, σ_true = spread of the simulated τ, ρ = σ_true / σ_fit(sim)",
                 "=" * 96]
        win = " | ".join(f"rms {lo:g}–{hi:g} s" for lo, hi in self.rms_windows).replace("–inf", "+")
        for t0 in self.compare_t_starts:
            S = self._res_stats.get(t0, [])
            if not S:
                continue
            g = lambda k: np.array([s_[k] for s_ in S], float)
            rms = np.array([s_["rms"] for s_ in S], float)
            lines.append(f"-- fit from {t0*1e3:.0f} ms --   [{win}]")
            if self.exp_y_full is not None:
                try:
                    R = self.fit_stats(x, self.exp_y_full, weights, t0)
                    lines.append(f"  REAL : τ = {R['tau']*1e3:8.3f} ms  σ_fit = {R['dtau']*1e3:.3f} ms  χ²ᵣ = {R['chi2r']:.2f}"
                                 f"  lag-1 = {R['lag1']:+.2f}  rms: " + "  ".join(f"{v:.2e}" for v in R["rms"]))
                except Exception as e:
                    R = None; lines.append(f"  REAL : fit failed ({e})")
            else:
                R = None
            lines.append(f"  SIM  : σ_true = {np.nanstd(g('tau'))*1e3:.3f} ms  σ_fit = {np.nanmean(g('dtau'))*1e3:.3f} ± "
                         f"{np.nanstd(g('dtau'))*1e3:.3f} ms  χ²ᵣ = {np.nanmean(g('chi2r')):.2f}  lag-1 = {np.nanmean(g('lag1')):+.2f}"
                         f"  rms: " + "  ".join(f"{m:.2e}" for m in np.nanmean(rms, 0)))
            if R is not None:
                ratio = R["dtau"] / np.nanmean(g("dtau"))
                lines.append(f"  σ_fit(real)/σ_fit(sim) = {ratio:.2f}   rms real/sim: "
                             + "  ".join(f"{a/b:.2f}" for a, b in zip(R["rms"], np.nanmean(rms, 0)))
                             + "   (ratio ≈ 1 = experiment correctly simulated)")
        lines.append("=" * 96)
        if self.boot_n > 0:
            lines += self.bootstrap_table(x, weights, tau_events_arr)
        txt = "\n".join(lines); print(txt)
        with open(self._outpath("comparison_real_vs_sim.txt"), "w") as fh:
            fh.write(txt + "\n")

    # ------------------------------------------------------------------
    # main loop
    # ------------------------------------------------------------------
    def start(self):
        logger.info(f'Simulation tasks: {self.settings_tasks}')

        # 'mle' reuses the actual decay times generated by 'addup' each trial (a genuine same-data benchmark),
        # so it cannot run on its own.
        if 'mle' in self.settings_tasks and 'addup' not in self.settings_tasks:
            logger.error("'mle' requires 'addup' to also be enabled (it reuses addup's per-trial decay times).")
            sys.exit(1)

        n_trials = self.params_n_trials
        # NaN, not 0: a failed trial must never enter the Gaussian fits as a fake data point
        tau_events_arr = np.full(n_trials, np.nan)
        sigma_events_arr = np.full(n_trials, np.nan)

        tau_distro_arr = np.full(n_trials, np.nan)
        sigma_distro_arr = np.full(n_trials, np.nan)

        tau_mle_arr = np.full(n_trials, np.nan)

        # bin starts 0, dt, 2dt, ... : exactly n_sim_steps bins
        x = np.arange(self.params_n_sim_steps) * self.params_timestep

        xvals = np.arange(
            -self.params_tau_seed, self.params_tau_seed, 2 * self.params_tau_seed / 1000
        )

        if 'addup' in self.settings_tasks:
            logger.info('Performing simulation task: addup')
            n_failed = 0
            weights = self.exp_s_full if self.fit_weights == "exp" else None
            for trial_number in tqdm(range(n_trials)):
                try:
                    y1, sigma_for_fit, trial_decay_times = self.create_from_events_with_fluctuations(
                        x, trial_number, add_empty_shots=self.params_empty_shots
                    )

                    # weights: the SAME for the real and the simulated fit (only their shape matters)
                    if weights is None and self.fit_weights == "model":
                        weights = sigma_for_fit                       # model sigma of the sum (first trial)
                    popt_events, pcov_events = self.fit_exponential(x, y1, weights)
                    # residual diagnostics for every fit start in compare_t_starts
                    for t0 in self.compare_t_starts:
                        self._res_stats.setdefault(t0, []).append(self.fit_stats(x, y1, weights, t0))
                    self._weights = weights
                    # bootstrap of the first bootstrap_trials simulated data sets (each one = one "experiment")
                    if self.boot_n > 0 and trial_number < self.boot_trials:
                        for t0 in self.compare_t_starts:
                            self._boot.setdefault(("sim", "res", t0), []).append(
                                self.bootstrap_residuals(x, y1, weights, t0))
                            self._boot.setdefault(("sim", "inj", t0), []).append(
                                self.bootstrap_injections(x, self._last_b_arr, weights, t0))
                    tau_events_arr[trial_number] = popt_events[1]
                    sigma_events_arr[trial_number] = np.sqrt(pcov_events[1, 1])  # error of tau from the fit

                    if 'mle' in self.settings_tasks:
                        tau_mle_arr[trial_number] = self.get_mle(trial_decay_times)

                    if self.settings_save_npz:
                        np.savez(self._outpath("exp_decay_sim.npz"), x=x, y1=y1, sigma=sigma_for_fit)

                    if self.settings_n_trial_plots is None or trial_number < self.settings_n_trial_plots:
                        self.plot_time_and_fit(
                            x,
                            y1,
                            popt_events,
                            limit=self.params_limit,
                            id_string=f"trial-{trial_number:04}_from_events",
                            display_fit=True,
                            pcov=pcov_events
                        )

                except (FloatingPointError, OptimizeWarning) as e:
                    logger.warning(e)
                    n_failed += 1
                    continue

            logger.info(f"{n_failed}/{n_trials} trials failed to fit and are excluded from the distributions.")
            if self.boot_n > 0:                                  # bootstrap of the REAL data
                logger.info(f"Bootstrap of the real data ({self.boot_n} resamples).")
                for t0 in self.compare_t_starts:
                    if self.exp_y_full is not None:
                        self._boot[("real", "res", t0)] = [self.bootstrap_residuals(x, self.exp_y_full, weights, t0)]
                if self.exp_injections_file:
                    x_i, Y_i = self.load_exp_injections()
                    w_i = weights if (weights is not None and len(weights) == len(x_i)) else None
                    for t0 in self.compare_t_starts:
                        self._boot[("real", "inj", t0)] = [self.bootstrap_injections(x_i, Y_i, w_i, t0)]
                        p_, c_ = self.fit_exponential(x_i, Y_i.sum(0), w_i, t0)
                        logger.info(f"  single-injection file, fit from {t0*1e3:.0f} ms: τ = {p_[1]*1e3:.3f} ± "
                                    f"{np.sqrt(c_[1,1])*1e3:.3f} ms (fit), injection bootstrap "
                                    f"{self._boot[('real','inj',t0)][0]*1e3:.3f} ms")
            self.print_comparison(x, weights, tau_events_arr)
            valid = np.isfinite(tau_events_arr) & np.isfinite(sigma_events_arr)
            if valid.sum() < 10:
                logger.error(f"Only {int(valid.sum())} valid trials -- too few for the Gaussian fits.")
                sys.exit(1)

            logger.info("Creating tau_tru distribution from add up spectra.")
            popt = self.fit_and_plot_gaussian(
                xvals, self.params_tau_seed - tau_events_arr[valid],
                id_string="tau_true_addup",
                label=r"$\tau_{\mathrm{true,addup}}$"
                #xvals, self.params_tau_seed - tau_events_arr[valid], id_string="$\tau_{\mathrm{true,addup}}$"
            )
            sigma_tru_from_addup = popt[2]
            logger.info(f"Sigma_tru from add up is = {sigma_tru_from_addup}")

            logger.info("Creating sigma distribution from add up spectra.")

            plt.figure()
            plt.hist(sigma_events_arr[valid], bins=50, alpha=0.7)
            plt.xlabel("σ_sim_old")
            plt.ylabel("Counts")
            plt.title("Distribution of σ_sim_old")
            plt.tight_layout()
            plt.savefig(self._outpath("sigma_sim_old_distribution.png"), dpi=150)
            plt.close()

            popt = self.fit_and_plot_gaussian(
                xvals, sigma_events_arr[valid],
                id_string="sigma_from_addup",
                label=r"$\sigma_{\mathrm{addup}}$"
                #xvals, sigma_events_arr[valid], id_string="$\sigma_{\mathrm{addup}}$"
            )
            mean_of_sigmas_from_addup = popt[1]
            logger.info(f"Mean of sigmas from addup is = {mean_of_sigmas_from_addup}")

            rho_from_addup = sigma_tru_from_addup / mean_of_sigmas_from_addup
            logger.info(f"rho from addup = {rho_from_addup}")

            logger.info("Creating the P distribion from add up.")
            self.fit_and_plot_gaussian(
                xvals, (self.params_tau_seed - tau_events_arr[valid]) / sigma_events_arr[valid],
                id_string="rho",
                label=r"$\rho$"
                #xvals, (self.params_tau_seed - tau_events_arr[valid]) / sigma_events_arr[valid], id_string="$\rho$"
            )

            if 'mle' in self.settings_tasks:
                logger.info("Creating tau_tru from unbinned MLE distribution (same decays as add-up).")
                valid_mle = np.isfinite(tau_mle_arr)
                popt_mle = self.fit_and_plot_gaussian(
                    xvals, self.params_tau_seed - tau_mle_arr[valid_mle],
                    id_string="tau_true_mle",
                    label=r"$\tau_{\mathrm{true,mle}}$"
                    #xvals, self.params_tau_seed - tau_mle_arr[valid_mle], id_string="$\tau_{\mathrm{true,mle}}$"
                )
                logger.info(f"Sigma_tru from MLE is = {popt_mle[2]}")

        if 'distro' in self.settings_tasks:
            logger.info('Performing simulation task: distro')
            for trial_number in tqdm(range(n_trials)):
                try:
                    y2 = self.create_from_distribution(x)
                    popt_distro, pcov_distro = self.fit_exponential(x, y2, None)
                    tau_distro_arr[trial_number] = popt_distro[1]
                    sigma_distro_arr[trial_number] = np.sqrt(pcov_distro[1, 1])  # error of tau from distro

                    self.plot_time_and_fit(
                        x, y2, popt_distro, self.params_limit, id_string=f'trial-{trial_number:04}_from_distro', display_fit=True, pcov=pcov_distro)

                except (FloatingPointError, OptimizeWarning) as e:
                    logger.warning(e)
                    continue

            valid_d = np.isfinite(tau_distro_arr) & np.isfinite(sigma_distro_arr)
            logger.info('Creating tau_tru from binned exponential distribution.')
            popt = self.fit_and_plot_gaussian(xvals, self.params_tau_seed - tau_distro_arr[valid_d], id_string='tau_tru_from_distro')
            sigma_tru_from_distro = popt[2]
            logger.info(f'Sigma_tru from distro is = {sigma_tru_from_distro}')

            logger.info('Creating sigma distribution from binned exponential distribution.')
            popt = self.fit_and_plot_gaussian(xvals, sigma_distro_arr[valid_d], id_string='sigma_from_distro')
            mean_of_sigmas_from_distro = popt[1]
            logger.info(f'Mean of sigmas from distro = {mean_of_sigmas_from_distro}')

            rho_from_distro = sigma_tru_from_distro / mean_of_sigmas_from_distro
            logger.info(f'rho from distro = {rho_from_distro}')

            logger.info('Creating the P distribion from distro.')
            self.fit_and_plot_gaussian(xvals, (self.params_tau_seed - tau_distro_arr[valid_d]) / sigma_distro_arr[valid_d], id_string='P_from_distro')

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
        outfilename = self._outpath(f"scatter_{self.simulation_duration:.2f}s")
        plt.savefig(outfilename + ".png", dpi=300)
        plt.savefig(outfilename + ".svg")
        plt.close()

        print("=" * 60)
        print("Aggregate plateau mean by occupancy (diagnostic check)")
        print("=" * 60)
        I0_check = self.params_mean_bkgnd
        A1_check = self.params_mean_ion - self.params_mean_bkgnd
        for n in sorted(self._diagnostic_level_by_n.keys()):
            values = self._diagnostic_level_by_n[n]
            observed_mean = np.mean(values)
            theoretical = I0_check + n * A1_check
            print(f"N={n}: observed avg={observed_mean:.4e}  theoretical={theoretical:.4e}  "
                  f"n_segments={len(values)}")
        #for late, rows in self._res_stats.items():
        #    r = np.array(rows); m, s = np.nanmean(r, 0), np.nanstd(r, 0)
        #    print(f"RionDS residuals, fit from {'167 ms' if late else 'first used bin'}: "
        #          f"rms 0–0.2 s {m[0]:.2e}±{s[0]:.1e} | 0.2–0.6 s {m[1]:.2e}±{s[1]:.1e} | 0.6+ s {m[2]:.2e}±{s[2]:.1e} | lag-1 {m[3]:+.2f}")


# -----------------

class Params(BaseModel):
    tau_seed: float
    timestep: float
    n_sim_steps: int
    n_sim: int
    n_trials: int
    n_decay_steps: int = Field(..., ge=0, le=1)
    mean_ion_number: float                      # non-integer means are allowed (Poisson mean)
    noise_model: str
    max_ions: int | None = None                 # optional: if set, clips the Poisson draws (lowers the realized mean)
    mean_ion: float
    stdv_ion: float
    mean_bkgnd: float
    stdv_bkgnd: float
    empty_shots: bool
    empty_shots_probability: float = Field(..., ge=0.0, le=1.0)
    limit: int
    seed: int | None = None
    ion_distribution: str = "poisson"          # "poisson" or "empirical"
    ion_distribution_file: str | None = None   # path to empirical .npy file containing initial ion numbers
    exp_decay_file: str | None = None
    fit_drop: list[int] = [0,1,2]              # bins excluded from the exponential fit (same as the real fit)
    fit_t_start: float = 0.0                   # [s] main fit uses only t >= fit_t_start (e.g. 0.167 as the notebook's 2nd start)
    fit_weights: str = "none"                  # "none" (no σ) | "model" (σ of the sum from the noise model) | "exp" (σ from exp_decay .npz)
    compare_t_starts: list[float] = [0.0, 0.167]   # [s] fit starts for the real-vs-sim comparison table
    rms_windows: list[list[float]] = [[0.0, 0.2], [0.2, 0.6], [0.6, 1e9]]   # [s] residual-rms windows (early/middle/tail)
    stdv_ion_corr: float = 0.0                 # optional time-correlated ion noise per ion (AR(1)); 0 = off
    corr_bins: float = 0.0                     # correlation length of that noise [bins]
    bootstrap_n: int = 0                       # bootstrap resamples (0 = off; 200 is enough for the σ to ~5 %)
    bootstrap_block: int = 1                   # block length [bins] of the residual bootstrap (1 = independent bins)
    bootstrap_trials: int = 10                 # simulated trials that are also bootstrapped (slow part)
    exp_injections_file: str | None = None     # npz with the REAL single-injection curves (Y, bins, dt) -> injection bootstrap
    rionsis_file: str | None = None            # RionSIS export (section 12): if set, ALL inputs and the real data come from it


class Settings(BaseModel):
    output_path: str
    plot_every_event: bool
    save_npz: bool
    tasks: list
    plot_titles: bool
    output_tag: str = ""                       # prefix of all output file names, e.g. "long_" / "short_"
    figure_formats: list[str] = ["png"]        # formats of the per-trial figures, e.g. ["pdf", "png"]
    n_trial_plots: int | None = None           # plot only the first N trials (None = all; plotting is the slow part)


class Config(BaseModel):
    params: Params
    settings: Settings


RIONSIS_KEYS = ("mean_ion_number", "mean_ion", "mean_bkgnd", "stdv_ion", "stdv_bkgnd")


def apply_rionsis_file(data: dict) -> dict:
    """rionsis_file set: ALL inputs come from that one RionSIS export (same injections, band, binning, units).
    Overrides the values in the toml (they are ignored) and logs every value taken from the file."""
    p = data.get("params", {})
    f = p.get("rionsis_file")
    if not f:
        return data
    d = np.load(f)
    take = {k: float(d[k]) for k in RIONSIS_KEYS}
    take["timestep"] = float(d["dt"])
    take["n_sim"] = int(d["n_inj"])
    take["n_sim_steps"] = int(d["Y"].shape[1])
    take["tau_seed"] = float(d["tau"])
    mu = take["mean_ion_number"]
    take["max_ions"] = int(np.ceil(mu + 6 * np.sqrt(mu)))
    take["exp_decay_file"] = f                 # summed curve + σ (keys y, sigma, first_bin)
    take["exp_injections_file"] = f            # single injections (keys Y, first_bin, dt)
    take["fit_drop"] = p.get("fit_drop", [])   # the first bins are already cut in RionSIS
    ignored = [k for k in take if k in p and k not in ("fit_drop",) and p[k] != take[k]]
    p.update(take)
    logger.info(f"rionsis_file {f}: all inputs taken from it ->")
    for k, v in take.items():
        if k not in ("exp_decay_file", "exp_injections_file"):
            logger.info(f"    {k} = {v}")
    if ignored:
        logger.info(f"    (values in the toml ignored for: {', '.join(ignored)})")
    data["params"] = p
    return data


def load_and_validate_toml(file_path: str) -> Config:
    data = apply_rionsis_file(toml.load(file_path))
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

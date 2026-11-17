#!/usr/bin/env python
"""
A particle decay time simulator for the investigation of the spectral add-up method

xaratustrah@github 2022-2026

Every injection: Poisson(mu) ions, exponential decay times (tau_seed), staircase I(n) = I0 + n*A1 with the
quadrature noise model sigma(n) = A1*sqrt(sigma0^2 + (k n)^2). The injections are summed and fitted with
A exp(-t/tau) + C, with the noise-model error bars (absolute_sigma=False, as the fit of the real data).
Repeated n_trials times: sigma_true = spread of tau, sigma_fit = fit error of tau, rho = sigma_true / sigma_fit.
"""

import argparse
import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import toml
from loguru import logger
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import ScaledTranslation
from pydantic import BaseModel, ConfigDict, Field, ValidationError
from scipy.optimize import OptimizeWarning, curve_fit
from tqdm import tqdm

# ======================================================================================================
# Figure format of the paper figures. FIG_WIDTH_IN = width in the LaTeX document (\the\linewidth / 72.27).
# ======================================================================================================
FIG_WIDTH_IN, FIG_HEIGHT_IN = 6.0, 3.6     # inches
FONT_AXIS, FONT_TICKS, FONT_LEGEND, FONT_ANNOT = 11, 10, 10, 8   # pt
LINEWIDTH, DPI = 1.5, 300
X_LEFT_S = -0.02                           # s, left x-limit of the time plots
GAMMA = 1.3957                             # Lorentz factor: lab -> rest frame
plt.rcParams.update({
    "font.size": FONT_AXIS, "axes.labelsize": FONT_AXIS, "axes.titlesize": FONT_AXIS,
    "xtick.labelsize": FONT_TICKS, "ytick.labelsize": FONT_TICKS, "legend.fontsize": FONT_LEGEND,
    "lines.linewidth": LINEWIDTH, "savefig.dpi": DPI, "pdf.fonttype": 42,
})
warnings.simplefilter("error", OptimizeWarning)
np.seterr(all="ignore")


# ======================================================================================================
# parameter file
# ======================================================================================================
class Params(BaseModel):
    model_config = ConfigDict(extra="forbid")
    # experiment (all from the RionSIS analysis, units of the RionSIS signal)
    tau_seed: float = Field(..., gt=0)          # [s] tau_lab of the fit of the real summed curve
    timestep: float = Field(..., gt=0)          # [s] time bin
    n_sim_steps: int = Field(..., gt=3)         # time bins
    n_sim: int = Field(..., gt=0)               # injections summed
    mean_ion_number: float = Field(..., gt=0)   # mu: mean ions per injection at t = 0 (Poisson)
    mean_bkgnd: float                           # I0: background level per injection (C / N)
    mean_ion: float                             # I1 = I0 + A1: level with one ion
    stdv_bkgnd: float = Field(..., ge=0)        # sigma(I0) = sigma0*A1: noise of an empty bin
    stdv_ion: float = Field(..., ge=0)          # sigma(I1) = k*A1: noise growth per ion
    ion_number_r: float = Field(0.0, ge=0)      # width r of a negative-binomial ion number (0 = Poisson)
    sigma_fit_real: float | None = None         # [s] fit error of tau of the real data (for rho and sigma_corr)
    fit_t_start: float = Field(0.0, ge=0)       # [s] first time used in the fit (first bin = 0), same as the real fit
    exp_data_file: str | None = None            # npz (t, y) of the real summed curve: overlaid in the plots
    # simulation
    n_decay_steps: int = Field(1, ge=0, le=1)   # 1: the bin of a decay holds the fraction the ion was alive
    empty_shots_probability: float = Field(0.0, ge=0, lt=1)   # extra empty injections (0 = none)
    n_trials: int = Field(1000, gt=9)
    seed: int | None = 1


class Settings(BaseModel):
    model_config = ConfigDict(extra="forbid")
    output_path: str = "./out"
    output_tag: str = ""                        # prefix of all output file names
    n_trial_plots: int = 3                      # summed curves plotted
    n_injection_plots: int = 10                 # single simulated injections plotted (first trial)
    plot_limit_injection: float | None = None   # [s] right end of the single-injection plots (None = whole record)
    plot_limit_sum: float | None = None         # [s] right end of the add-up plots (None = whole record)
    figure_formats: list[str] = ["png"]
    save_npz: bool = False


class Config(BaseModel):
    model_config = ConfigDict(extra="forbid")
    params: Params
    settings: Settings = Settings()


# ======================================================================================================
class Simulator:
    def __init__(self, cfg: Config):
        p, s = cfg.params, cfg.settings
        self.p, self.s = p, s
        self.tau = p.tau_seed
        self.dt = p.timestep
        self.N = p.n_sim
        self.I0, self.A1 = p.mean_bkgnd, p.mean_ion - p.mean_bkgnd     # A1 = I1 - I0: single-ion step
        self.s_bkg, self.s_ion = p.stdv_bkgnd, p.stdv_ion
        self.x = np.arange(p.n_sim_steps) * self.dt  # bin starts
        self.limit_s = self.x[-1] + self.dt
        self.y_exp = self.fit_exp = None
        if p.exp_data_file:
            d = np.load(p.exp_data_file)
            if len(d["y"]) != p.n_sim_steps:
                logger.warning(f"exp_data_file has {len(d['y'])} bins, n_sim_steps = {p.n_sim_steps}: not overlaid")
            else:
                self.y_exp = np.asarray(d["y"], float)
        self.rng = np.random.default_rng(p.seed)
        self.out = os.path.join(s.output_path, "")
        os.makedirs(self.out, exist_ok=True)
        # error bars of the summed curve: noise model with the expected occupancy (same as RionSIS):
        # σ_sum(t)² = N·(σ(I0)² + σ(I1)²·E[n(t)²]),  E[n²] = (1-p_empty)·(μ(t) + μ(t)²),  μ(t) = μ·exp(-t/τ)
        # (= N·A1²·(σ0² + k²·E[n²]) with σ0 = σ(I0)/A1, k = σ(I1)/A1)
        mu_t = p.mean_ion_number * np.exp(-self.x / self.tau)
        n2 = (1 - p.empty_shots_probability) * (mu_t + mu_t**2 * (1 + (1 / p.ion_number_r if p.ion_number_r else 0)))
        self.sigma_sum = np.sqrt(self.N * (self.s_bkg**2 + self.s_ion**2 * n2))
        logger.info(f"output: {self.out}{s.output_tag}*   seed: {p.seed}")

    # ---------------------------------------------------------------- helpers
    @staticmethod
    def expo(t, A, tau, C):
        return A * np.exp(-t / tau) + C

    @staticmethod
    def gauss(x, a, m, s):
        return a * np.exp(-((x - m) ** 2) / (2.0 * s**2))

    def path(self, stem):
        return f"{self.out}{self.s.output_tag}{stem}"

    def save(self, fig, stem):
        for ext in self.s.figure_formats:
            fig.savefig(f"{self.path(stem)}.{ext}", dpi=DPI)
        plt.close(fig)

    def style(self, ax, limit=None):
        ax.set_xlabel("Time [s]"); ax.set_ylabel("Amplitude [a.u.]")
        ax.set_xlim(X_LEFT_S, limit or self.limit_s)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(FONT_TICKS)
        ax.grid(True)
        ax.legend(fontsize=FONT_LEGEND, loc="upper right")

    def fit(self, y):
        """A exp(-t/tau) + C with the noise-model error bars, errors scaled by sqrt(chi2_r)."""
        m = self.x >= self.p.fit_t_start - 1e-9
        x, y_, s = self.x[m], y[m], self.sigma_sum[m]
        p0 = [max(y_[0] - y_[-1], 1e-15) * np.exp(x[0] / self.tau), self.tau, y_[-1]]
        popt, pcov = curve_fit(self.expo, x, y_, p0=p0, sigma=s, absolute_sigma=False,
                               bounds=([0, 1e-4, -np.inf], [np.inf, 50, np.inf]), maxfev=20000)
        r = (y_ - self.expo(x, *popt)) / s
        return popt, pcov, np.sum(r**2) / (len(r) - 3)

    # ---------------------------------------------------------------- one trial
    def simulate_trial(self, trial):
        """n_sim injections, summed. Returns (sum, all decay times)."""
        x, n_bins = self.x, len(self.x)
        y_sum = np.zeros(n_bins)
        decays = []
        for i in range(self.N):
            r_nb = self.p.ion_number_r
            n_ions = int(self.rng.poisson(self.rng.gamma(r_nb, self.p.mean_ion_number / r_nb) if r_nb
                                          else self.p.mean_ion_number))
            if self.p.empty_shots_probability > 0 and self.rng.random() < self.p.empty_shots_probability:
                n_ions = 0                                # forced empty shot
            if n_ions == 0:
                b = self.rng.normal(self.I0, self.s_bkg, n_bins)
                n_alive, t_dec = np.zeros(n_bins, int), np.array([])
            else:
                t_dec = np.sort(self.rng.exponential(self.tau, n_ions))
                decays.append(t_dec)
                n_alive = n_ions - np.searchsorted(t_dec, x, side="right")
                level = self.I0 + n_alive * self.A1
                if self.p.n_decay_steps == 1:            # decay bin: ion contributes only while alive
                    for td in t_dec:
                        j = int(td / self.dt)
                        if j < n_bins:
                            level[j] -= (1 - (td % self.dt) / self.dt) * self.A1
                sigma = np.sqrt(self.s_bkg**2 + (n_alive * self.s_ion) ** 2)
                b = level + self.rng.normal(0, sigma)
            y_sum += b
            if trial == 0 and i < self.s.n_injection_plots:
                self.plot_injection(n_alive, b, i)
        return y_sum, (np.concatenate(decays) if decays else np.array([]))

    def mle(self, t_dec):
        """Censored MLE of tau from the same decay times (best possible precision)."""
        T = self.limit_s
        seen = t_dec[t_dec < T]
        return np.nan if len(seen) == 0 else seen.mean() + (len(t_dec) - len(seen)) * T / len(seen)

    # ---------------------------------------------------------------- figures
    def plot_injection(self, n_alive, b, i):
        lim = self.s.plot_limit_injection or self.limit_s
        keep = self.x < lim                                   # cut the trace at the plot limit
        x, b, n_alive = self.x[keep], b[keep], n_alive[keep]
        fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), constrained_layout=True)
        ax.stairs(b, np.append(x, x[-1] + self.dt), color="C0", label="Simulated signal")   # one step per time bin
        ch = np.where(np.diff(n_alive) != 0)[0]
        bounds = np.concatenate(([0], ch + 1, [len(n_alive)]))
        for j in range(len(bounds) - 1):
            a, e = bounds[j], bounds[j + 1]
            n = n_alive[a]
            lev = self.I0 + n * self.A1
            ax.hlines(lev, x[a], min(x[e - 1] + self.dt, lim), colors="red", linestyle="--", alpha=0.6,
                      label="Staircase level" if j == 0 else None)
            ax.text((x[a] + min(x[e - 1] + self.dt, lim)) / 2, lev, f"{n} ion{'s' if n != 1 else ''}",
                    transform=ax.transData + ScaledTranslation(0, 3 / 72, fig.dpi_scale_trans), ha="center",
                    va="bottom", fontsize=FONT_ANNOT, bbox=dict(facecolor="white", alpha=0.8, pad=1.0))
        for t in x[ch]:
            ax.axvline(t, color="grey", linestyle="--", linewidth=1)
        self.style(ax, lim)
        self.save(fig, f"multi_decay_trial_{i:04}")

    def plot_sum(self, y, popt, pcov, trial):
        lim = self.s.plot_limit_sum or self.s.plot_limit_injection or self.limit_s   # same cut as the staircases
        k = self.x < lim
        x = self.x[k]
        fig, (ax, axr) = plt.subplots(2, 1, sharex=True, figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN * 1.3),
                                      gridspec_kw=dict(height_ratios=[3, 1]), constrained_layout=True)
        edges = np.append(x, x[-1] + self.dt)
        ax.stairs(y[k], edges, color="C0", label="Summed-up simulated signal")
        axr.stairs(y[k] - self.expo(x, *popt), edges, color="C0")                 # simulated - simulated fit
        if self.y_exp is not None:
            ax.stairs(self.y_exp[k], edges, color="C1", alpha=0.6, label="Experimental data")
            if self.fit_exp is None:
                self.fit_exp = self.fit(self.y_exp)[0]                             # same fit on the real curve
            axr.stairs(self.y_exp[k] - self.expo(x, *self.fit_exp), edges, color="C1", alpha=0.6)
        tf = np.linspace(self.p.fit_t_start, x[-1], 400)
        ax.plot(tf, self.expo(tf, *popt), "--", color="red", alpha=0.6, label="Simulated fit")
        if self.p.fit_t_start > 0:
            for a in (ax, axr):
                a.axvline(self.p.fit_t_start, color="k", ls=":", lw=1)
        self.style(ax, lim)
        ax.set_xlabel("")
        ax.yaxis.set_major_locator(MultipleLocator(10 ** np.floor(np.log10(np.max(np.abs(y)))) / 2))
        axr.axhline(0, color="k", lw=0.8)
        axr.set_xlabel("Time [s]"); axr.set_ylabel("Residuals", fontsize=FONT_LEGEND)
        axr.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axr.yaxis.get_offset_text().set_fontsize(FONT_TICKS)
        axr.grid(True)
        self.save(fig, f"trial-{trial:04}_summed")
        if self.s.save_npz:
            np.savez(self.path(f"trial-{trial:04}_summed.npz"), x=x, y=y, sigma=self.sigma_sum,
                     yfit=self.expo(x, *popt))

    def plot_hist(self, vals, stem, label, xlabel):
        """Histogram with the Gaussian of its mean and standard deviation (the same numbers as in summary.txt)."""
        vals = np.asarray(vals)[np.isfinite(vals)]
        cnt, edges = np.histogram(vals, bins=100)
        xc = (edges[:-1] + edges[1:]) / 2
        m, s = float(np.mean(vals)), float(np.std(vals, ddof=1))
        a = len(vals) * (edges[1] - edges[0]) / (s * np.sqrt(2 * np.pi))       # same area as the histogram
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(vals, bins=100, alpha=0.6, label=label)
        xf = np.linspace(xc[0], xc[-1], 400)
        ax.plot(xf, self.gauss(xf, a, m, s), "r", lw=2,
                label="Gaussian\n" + fr"$\mu =$ {m:0.2e}" + "\n" + fr"$\sigma =$ {s:.2e}")
        ax.axvline(m, color="k", linestyle="--", linewidth=1.2)
        ax.annotate("", xy=(m + s, a * np.exp(-0.5)), xytext=(m, a * np.exp(-0.5)),
                    arrowprops=dict(arrowstyle="<->", color="k", lw=1.5))
        ax.text(m + s / 2, a * np.exp(-0.5), r"$\sigma$", ha="center", va="bottom", fontsize=16)
        ax.set_xlabel(xlabel); ax.set_ylabel("Counts"); ax.legend(fontsize=16, loc="upper right")
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(f"{self.path(stem)}.{ext}", dpi=300)
        plt.close(fig)
        return m, s

    # ---------------------------------------------------------------- main loop
    def run(self):
        nT = self.p.n_trials
        tau, dtau, chi2, tau_mle = (np.full(nT, np.nan) for _ in range(4))
        for t in tqdm(range(nT)):
            y, t_dec = self.simulate_trial(t)
            tau_mle[t] = self.mle(t_dec)
            try:
                popt, pcov, chi2[t] = self.fit(y)
            except (RuntimeError, OptimizeWarning, ValueError) as e:
                logger.warning(f"trial {t}: fit failed ({e})")
                continue
            tau[t], dtau[t] = popt[1], np.sqrt(pcov[1, 1])
            if t < self.s.n_trial_plots:
                self.plot_sum(y, popt, pcov, t)
        ok = np.isfinite(tau)
        if ok.sum() < 10:
            logger.error(f"only {ok.sum()} valid trials"); sys.exit(1)
        logger.info(f"{nT - ok.sum()}/{nT} trials failed and are excluded.")

        # distributions (paper figures)
        _, sig_gauss = self.plot_hist(self.tau - tau[ok], "tau_true_addup", r"$\tau_{\mathrm{true,addup}}$",
                                      r"$\tau_{seed} - \tau$ [s]")
        self.plot_hist(dtau[ok], "sigma_from_addup", r"$\sigma_{\mathrm{addup}}$", r"$\sigma_{\mathrm{fit}}$ [s]")
        self.plot_hist((self.tau - tau[ok]) / dtau[ok], "pull", "pull", r"$(\tau_{seed} - \tau)/\sigma_{\mathrm{fit}}$")
        # ρ is one number per run: its distribution comes from resampling the trials (Monte-Carlo uncertainty of ρ)
        t_ok, d_ok = tau[ok], dtau[ok]
        idx = self.rng.integers(0, len(t_ok), (len(t_ok), len(t_ok)))
        rho_rs = np.std(t_ok[idx], axis=1, ddof=1) / np.mean(d_ok[idx], axis=1)
        self.plot_hist(rho_rs, "rho", r"$\rho$ (resampled trials)", r"$\rho = \sigma_{\mathrm{true,addup}}/\sigma_{\mathrm{addup}}$")
        self.rho_err = float(np.std(rho_rs, ddof=1))
        _, sig_mle_g = self.plot_hist(self.tau - tau_mle[np.isfinite(tau_mle)], "tau_true_mle",
                                      r"$\tau_{\mathrm{true,mle}}$", r"$\tau_{seed} - \tau_{\mathrm{MLE}}$ [s]")
        self.summary(tau[ok], dtau[ok], chi2[ok], tau_mle[np.isfinite(tau_mle)], sig_gauss)

    def summary(self, tau, dtau, chi2, tau_mle, sig_gauss):
        p, ms = self.p, 1e3
        s_true, s_fit, s_fit_sd = np.std(tau, ddof=1), np.mean(dtau), np.std(dtau, ddof=1)
        rho = s_true / s_fit
        hl = np.log(2) / GAMMA
        L = ["=" * 92,
             f"RionDS  {p.n_trials} trials | N = {p.n_sim}, μ = {p.mean_ion_number:.4f}, I0 = {self.I0:.4e}, "
             f"I1 = {p.mean_ion:.4e}, σ(I0) = {self.s_bkg:.4e} ({self.s_bkg/self.A1:.4f} A1), "
             f"σ(I1) = {self.s_ion:.4e} ({self.s_ion/self.A1:.4f} A1)\n  intermediate frame = {bool(p.n_decay_steps)}, "
             f"empty shots = {p.empty_shots_probability:g}, ion number "
             f"{'negative binomial r = %g' % p.ion_number_r if p.ion_number_r else 'Poisson'}",
             f"  τ_seed = {self.tau*ms:.3f} ms, {p.n_sim_steps} bins × {self.dt*ms:.3f} ms, fit with noise-model σ "
             f"from {p.fit_t_start*ms:.1f} ms",
             "=" * 92,
             f"  τ_lab (sim, mean)             {np.mean(tau)*ms:9.3f} ms   (bias {np.mean(tau-self.tau)*ms:+.3f} ms)",
             f"  σ_true  (std of τ)            {s_true*ms:9.3f} ms",
             f"  σ_fit   (mean fit error)      {s_fit*ms:9.3f} ± {s_fit_sd*ms:.3f} ms   χ²ᵣ = {np.mean(chi2):.2f}",
             f"  σ_MLE   (unbinned, same decays){np.std(tau_mle, ddof=1)*ms:8.3f} ms",
             f"  ρ = σ_true / σ_fit            {rho:9.3f} ± {getattr(self, 'rho_err', np.nan):.3f} (MC)"]
        if p.sigma_fit_real is not None:
            sr = p.sigma_fit_real
            L += ["-" * 92,
                  f"  REAL σ_fit = {sr*ms:.3f} ms  ->  real / sim = {sr/s_fit:.2f}   "
                  f"(pull {(sr-s_fit)/s_fit_sd:+.1f} σ; |pull| < 2: experiment correctly simulated)",
                  f"  σ_corr = σ_fit(real) × ρ = {sr*rho*ms:.3f} ms",
                  f"  τ_lab = {self.tau*ms:.3f} ± {sr*rho*ms:.3f} ms   ->   T½ (rest) = "
                  f"{self.tau*hl*ms:.3f} ± {sr*rho*hl*ms:.3f} ms"]
        L.append("=" * 92)
        txt = "\n".join(L)
        print(txt)
        with open(self.path("summary.txt"), "w") as fh:
            fh.write(txt + "\n")


# ======================================================================================================
def main():
    logger.remove(0)
    logger.add(sys.stdout, level="INFO")
    ap = argparse.ArgumentParser(prog="decaysim")
    ap.add_argument("param", help="parameter file (.toml)")
    args = ap.parse_args()
    try:
        cfg = Config(**toml.load(args.param))
    except ValidationError as e:
        logger.error(f"parameter file {args.param}:\n{e}")
        sys.exit(1)
    try:
        Simulator(cfg).run()
    except KeyboardInterrupt:
        logger.info("aborted by user")


if __name__ == "__main__":
    main()

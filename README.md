# decaysim

A decay simulator.


## Motivation


## Installation

You can directly use the app:

```
python3 -m decaysim parameters.toml
```

or first install it by:

```
pip install -r requirements.txt
pip install .
```

Uninstalling is similarly done by calling `pip`.

## Simulation settings

The simulation is setup using a toml file. The following are the description of the parameters:

| Parameter | Type | Range | Example Value | Description |
|-----------|------|-------|---------------|-------------|
| tau_seed | float | Any positive number | 23.08 | Estimated decay lifetime |
| timestep | float | Any positive number | 0.083886 | Length of one time frame within which a decay can happen |
| n_sim_steps | int | Any positive integer | 477 | Observation duration, usually an integer multiple of timestep. It is suggested that this time is long enough for most of the particles to decay.|
| n_sim | int | Any positive integer | 32 | Number of decays to be simulated or injections |
| n_trials | int | Any positive integer | 1 | Experiment repetition count (dice throws) |
| n_decay_steps | int | 0 or 1 | 0 | Indicates 0 (i.e. no) intermediate step, or 1 step in the decay |
| mean_ion | float | Any positive number | 9.44e-7 | Amplitude of a single ion |
| stdv_ion | float | Any positive number | 7.7408e-8 | Error in ion amplitude |
| mean_bkgnd | float | Any positive number | 5.44e-7 | Mean background noise level |
| stdv_bkgnd | float | Any positive number | 2.06e-8 | Error in background noise level |
| empty_shots | bool | true or false | false | Includes empty shots in simulation |
| empty_shots_probability | float | 0.0 - 1.0 | 0.0 | Probability of an empty shot occurring |
| output_path | str | Valid directory path | "." | Directory for storing results |
| plot_every_event | bool | true or false | false | Determines if every event should be plotted (high computation cost) |


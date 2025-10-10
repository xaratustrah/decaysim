# decaysim

A particle decay time simulator for the investigation of the spectral add up method


## Motivation

Due to their non-destructive nature, Schottky detectors have been preferably used in storage ring experiments in order to determine lifetimes of heavy ion beams in a variety of experiments. With the development of more sensitive detectors, it has been possible to increase the detection sensitivity down to a few or even single ions. This is crucial in experiments involving low-yield nuclei. The increased sensitivity also allows for faster detection, thanks to the reduced recording and/or averaging time needed for a sufficient signal to noise ratio (SNR). This is useful for lifetime measurements of short-lived nuclear species and their isomeric states. Targeting these extreme experimental scenarios, this work presents the spectral add-up technique, which, when combined with the isochronous ion optical mode of the heavy ion Experimental Storage Ring (ESR), can significantly increase the accuracy of lifetime measurements. The efficacy of the method is studied using dedicated simulations.

For more information please refer to the publication related to this software.

## Installation and usage

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

#### Simulation settings

The simulation is setup using a [TOML](https://toml.io/en/) configuration file. The following are the description of the parameters:

| Parameter | Type | Range | Example Value | Description |
|-----------|------|-------|---------------|-------------|
| tau_seed | float | Any positive number | 33 | Estimated decay lifetime |
| timestep | float | Any positive number | 0.671 | Length of one time frame within which a decay can happen |
| n_sim_steps | int | Any positive integer | 210 | Observation duration, usually an integer multiple of timestep. It is suggested that this time is long enough for most of the particles to decay.|
| n_sim | int | Any positive integer | 32 | Number of decays to be simulated or injections |
| n_trials | int | Any positive integer | 1000 | Experiment repetition count (dice throws) |
| n_decay_steps | int | 0 or 1 | 0 | Indicates 0 (i.e. no) intermediate step, or 1 step in the decay |
| mean_ion | float | Any positive number | 9.44e-7 | Amplitude of a single ion |
| stdv_ion | float | Any positive number | 7.7408e-8 | Error in ion amplitude |
| mean_bkgnd | float | Any positive number | 5.44e-7 | Mean background noise level |
| stdv_bkgnd | float | Any positive number | 2.06e-8 | Error in background noise level |
| empty_shots | bool | true or false | false | Includes empty shots in simulation |
| empty_shots_probability | float | 0.0 - 1.0 | 0.0 | Probability of an empty shot occurring |
| tasks | list | 'addup', 'distro', 'mle' | 'addup' | Which method to simulate |
| output_path | str | Valid directory path | "./out" | Directory for storing results |
| plot_every_event | bool | true or false | false | Determines if every event should be plotted (high computation cost) |
| save_npz | bool | true or false | false | Whether or not to save data as binary files |
| plot_titles | bool | true or false | false | Whether plots should have a title |



#### Automatic rerun of the simulation

Although many measures have been implemented to prevent crashes, the simulation can sometimes fail due to incorrect parameters. You can simply delete the old results and start again. It is also possible to automate this process using a bash script like the following: This assumes that the output directory is set to 'out' in the parameters file.

```bash
#!/bin/bash
while ! python3 -m decaysim simparams_72br.toml; do
    echo "Oops, the simnulation didn't work!"
    rm ./out/*
done
echo "Simulation finished at $(date)" > run.log

```

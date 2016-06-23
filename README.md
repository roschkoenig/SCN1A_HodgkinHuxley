## Modelling dynamic abnormalities caused by an SCN1A mutation
The Matlab scripts uploaded here implement a Hodgkin-Huxley type neuronal model with parameterisation derived directly from experimental voltage clamp results. They correspond in part to the <a href="http://dx.doi.org/10.1101/048520">preprint</a> of an upcoming publication. The scripts were used to produce the figures in the final manuscript, which will be published soon:

Peters C<sup>+</sup>, Rosch RE<sup>+</sup>, Hughes E, Ruben PC (2016): Temperature-dependent changes in neuronal dynamics in a patient with an SCN1A mutation and hyperthermia induced seizures. <em> under review </em> <br>
<sup>+</sup> <em>equal contribution</em>

### Fitting Boltzmann Equations and Plotting Steady-State Parameters
```
steady_state_curves
```
![Fitting Steady State Parameters](https://cloud.githubusercontent.com/assets/12950773/16322255/bdec9b9a-3999-11e6-8d54-b0e5a8662018.png)

In the first instance, Boltzmann equations were fitted to an existing formulation of the Hodgkin Huxley Model that describes the dynamics at the membrane of a mammalian cortical neuron. The first component of the script `steady_state_curves` fits the Boltzmann equations. 
This formulation is taken as the baseline for all further modelling - i.e. the wildtype at 37ºC is modelled as this baseline parameter composition. <br>

![Steady State Parameters for Different Conditions](https://cloud.githubusercontent.com/assets/12950773/16322251/b98c8056-3999-11e6-939d-44f6979896f8.png)

Remaining empirical values are represented as absolute (voltage parameters), or relative (slope parameters) deviations from the baseline in the model. The steady states that this produces for the sodium channel gating parameters are plotted in the second half of the `steady_state_curves` script. 

### Virtual Voltage Clamp to Assess for Fast Inactivation
``` 
cort_fastinactivation.m
```
![Fast Inactivation](https://cloud.githubusercontent.com/assets/12950773/16312374/505c5d86-396c-11e6-98c9-e8e1f5e77388.png)

This function evaluates model dynamics by performing a 'virtual voltage clamp' experiments for the four parameterisations of the model (i.e. baseline, wild type sodium channels at 40ºC, SCN1A mutation at 37 and 40ºC). This consists of voltage clamp (i.e. specific voltages being enforced on the model), delivering paired voltage pulses with increasing inter-pulse intervals. Changes in sodium conductance are estimated in a reduced Hodgkin Huxley Model implemented in `cort_hh_fastinactivation.m`. Current peaks resulting from these paired stimuli are estimated and plotted for each of the individual parameterisations of the model in the figure. 


### Running Models Under Different Conditions
```
cort_variable_hh_eval
```
![Output from full model simulations](https://cloud.githubusercontent.com/assets/12950773/16312193/8040b0fc-396b-11e6-8560-272a28e32430.png "Model simulations")
This script calls `cort_variable_hh.m` and uses Matlab ODE solver to integrate the model and plot the time course. The parameters are derived from a baseline model of cortical dynamics with experimental values implemented in terms of their deviations from those baseline parameters. 

### Evaluation of Model Bifurcations 
```
cort_bifurcation.m
```

![Bifurcation Analysis](https://cloud.githubusercontent.com/assets/12950773/16322259/c2c101e2-3999-11e6-9547-5cf1b7d2247d.png)
This routines calls `cort_variable_hh.m` and to run simulations of the cortical Hodgkin Huxley model across the four different parameterisations (i.e. baseline, wild type sodium channels at 40ºC, SCN1A mutation at 37 and 40ºC). The model is estimated to steady state, and final state values from one simulation are used as initial values for the next simulation. This simulation is performed repeatedly across increasing, and then decreasing values of the input current Parameter I<sub>stim</sub>. The steady state values of each model run are then plotted as bifurcation plots, separately for the four parameterisations of the model.

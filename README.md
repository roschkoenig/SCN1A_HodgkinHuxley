## Modelling dynamic abnormalities caused by an SCN1A mutation
The Matlab scripts uploaded here implement a Hodgkin-Huxley type neuronal model with parameterisation derived directly from experimental voltage clamp results. They correspond in part to the <a href="http://dx.doi.org/10.1101/048520">preprint</a> of an upcoming publication. The scripts were used to produce the figures in the final manuscript, which will be published soon:

Peters C<sup>+</sup>, Rosch RE<sup>+</sup>, Hughes E, Ruben PC (2016): Temperature-dependent changes in neuronal dynamics in a patient with an SCN1A mutation and hyperthermia induced seizures. <em> under review </em> <br>
<sup>+</sup> <em>equal contribution</em>

### Fitting Boltzmann Equations and Plotting Steady-State Parameters
```
steady_state_curves
```
In the first instance, Boltzmann equations were fitted to an existing formulation of the Hodgkin Huxley Model that describes the dynamics at the membrane of a mammalian cortical neuron. The first component of the script `steady_state_curves` fits the Boltzmann equations. 
This formulation is taken as the baseline for all further modelling - i.e. the wildtype at 37ºC is modelled as this baseline parameter composition. <br>
Remaining empirical values are represented as absolute (voltage parameters), or relative (slope parameters) deviations from the baseline in the model. The steady states that this produces for the sodium channel gating parameters are plotted in the second half of the `steady_state_curves` script. 

### Running Models Under Different Conditions
```
cort_variable_hh_eval
```
![Output from full model simulations](https://cloud.githubusercontent.com/assets/12950773/16312193/8040b0fc-396b-11e6-8560-272a28e32430.png "Logo Title Text 1")
This script calls `cort_variable_hh.m` and uses Matlab ODE solver to integrate the model and plot the time course. The parameters are derived from a baseline model of cortical dynamics with experimental values implemented in terms of deviations from the baseline parameters. 



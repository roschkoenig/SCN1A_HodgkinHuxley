%% Conducting a virtua experiment estimating fast inactivation
%==========================================================================
% The appropriateness of the model is evaluated by conducting virtual
% voltage clamp experiments in order to evaluate whether the model predicts
% real empirical measures. 
% The routine below performs an two-pulse voltage clamp experiment aimed to
% evaluate the dynamics of fast inactivation. A reduced set of HH equations
% only including sodium channel function is probed with a specific set of
% paired input pulses with increasing inter-pulse intervals. Short
% intervals leave little time for inactivation to recover, and therefore
% diminish model output.
% The model performs the experiment across a range of inter-pulse intervals
% and plots the maximally achieved pulse amplitude of the second pulse for
% a wild type sodium channel (baseline), as well as a specific SCN1A
% mutation causing epilepsy. 
%
% The routine produces a figure that was the basis for Figure 5b and d in 
% the publication below:
% 
% Peters C, Rosch RE, Hughes E, Ruben P (2016) Temperature-dependent 
% changes in neuronal dynamics in a patient with an SCN1A mutation and 
% hyperthermia induced seizures (under review)

% Housekeeping
%==========================================================================
clear all
clf
count = 1;

% Manual parameters
%==========================================================================
% Define manually the name of the experiments and which experimental
% parameters are going to be modelled in this evaluation

% Choose experiments to be modelled
%--------------------------------------------------------------------------
n_exp           = {'WT37','WT40','AV37','AV40'};
exp_parameters  = [1 2 3 4];        % 1 = WT37 (norm), 2 = WT40, 3 = AV37, 3 = AV40;
plot_phase      = 0;
I_stim          = 0;                % I_stim    HH value = 0.200 nA

% Set parameters (currently done according to literature
%--------------------------------------------------------------------------
params(1) = 0.010;      % C         HH
params(2) = 0.0000205; 	% g_L       Pospischil 2008, Fig 2
params(3) = 5;          % g_K       Traub 1991
params(4) = 56;         % g_Na      Pospischil 2008, Fig 2
params(5) = -90;        % E_K       Traub 1991
params(6) = -70.3;    	% E_L       Pospischil 2008, Fig 2
params(7) = 50;         % E_Na      Traub 1991
params(11) = -60;      	% V_t       arbitrary (guided by Posposchil)
params(8) = I_stim;  

for e = 1:4; %exp_parameters
  
% Parameters defined according to experimental values (with WT37 as baseline)
%==========================================================================
switch e
    case 1  % standard HH formulation (WT37)
    %----------------------------------------------------------------------
        temperature = 37 + 273;       
        dV2_m   = 0;        % normalised to 0
        dV2_h   = 0;        % normalised to 0
        t_off   = 0;        % normalised to 0
        m_z     =  3.6089;  % normalised to cortical neuron values
        h_z     = -6.6210;  % normalised to cortical neuron values
        
    case 2  % Experimental parameters for WT40
    %----------------------------------------------------------------------
        temperature = 40 + 273;
        dV2_m   =  6.49;    % absolute difference to WT37
        dV2_h   = -1.9;     % absolute difference to WT37
        t_off   = 1.9;      % absolute difference to WT37
        m_z    	=  2.43*3.6089/2.14;    % relative difference to WT37
        h_z   	= -3.79*6.6210/3.14;    % relative difference to WT37
        
    case 3  % Experimental parameters for AV37
    %----------------------------------------------------------------------
        temperature = 37 + 273;
        dV2_m   =  5.94;    % absolute difference to WT37
        dV2_h   =  2.4;     % absolute difference to WT37
        t_off   = -2.4;     % absolute difference to WT37
        m_z    	=  2.58*3.6089/2.14;    % relative difference to WT37     
        h_z   	= -3.73*6.6210/3.14;    % relative difference to WT37
        
    case 4  % Experimental parameters for AV40
    %----------------------------------------------------------------------
        temperature = 40 + 273;
        dV2_m   = 9.93;     % absolute difference to WT37
        dV2_h   = 11.7;     % absolute difference to WT37
        t_off   = -11.7;    % absolute difference to WT37
        m_z    	= 2.1*3.6089/2.14;  	% relative difference to WT37 
        h_z   	= -3.38*6.6210/3.14;    % relative difference to WT37    
end

% Calculate model parameters from values given above
%--------------------------------------------------------------------------
s_m         = (0.0863 * temperature)/m_z;   % slope derived from m_z
s_h         = -(0.0863 * temperature)/h_z;  % slope derived from h_z

params(9)   = -39.0 + dV2_m;	% difference from standard cortical neuron
params(10)  = s_m;           
params(12)  = -43.2 + dV2_h;    % difference from standard cortical neuron
params(13)  = s_h;    

params(14) = t_off;         


g_Na        = params(4);
E_Na        = params(7);
x_ini       = [0 0 0 0];

% Run Model
%==========================================================================
rec_t(1) = 0.1;      % recovery time in miliseconds
for r = 1:10
p1          = [100 200];
p2          = [0 10] + [p1(2) + rec_t(r)];
t_range     = [0 p2(2)+100];

options     = odeset('InitialStep',0.005,'MaxStep',t_range(2)/1000);
[t,x]       = ode45(@(t,x)cort_hh_fastinactivation(t,x,params, p1, p2),t_range,x_ini,options);

% Estimate currents for each time point
%--------------------------------------------------------------------------
% Estimate actual currents from model output (sodium channel conductance) 
% and the set voltages as per voltage clamp experimental design

for s = 1:length(t)
    
    if t(s) < p1(1), vol(s) = -130;                      % Initial cond
    elseif t(s) > p1(1) && t(s) < p1(2), vol(s) = 0;     % Pulse 1
    elseif t(s) > p2(1) && t(s) < p2(2), vol(s) = 0;     % Pulse 2 
    else vol(s) = -90;                                   % Recovery Time
    end
    
    x(s,3) = (-g_Na * x(s,1)^3 * x(s,2) * (vol(s) - E_Na));    
end

% Estimate current peak after second pulse
%--------------------------------------------------------------------------
count = 1;
for i = 1:length(t)
    if t(i) > p2(1) && t(i) < (p2(1) + 10)
        ind(count) = i;
        count = count+1;
    end
end

[val loc]    = findpeaks(x(ind,3));
IR(e,r,1)    = rec_t(r);
IR(e,r,2)    = max(val);

rec_t(r+1)   = rec_t(r) * 2;
end
clear rec_t p1 p2 ind
end

subplot(2,1,1)
plot(log(IR(1,:,1)), IR(1,:,2)/max(IR(1,:,2)), 'ks-'); hold on
plot(log(IR(3,:,1)), IR(3,:,2)/max(IR(3,:,2)), 'ks:')
title('37 degrees');
legend({'WT', 'AV'});

subplot(2,1,2)
plot(log(IR(2,:,1)), IR(2,:,2)/max(IR(2,:,2)), 'rs-'); hold on
plot(log(IR(4,:,1)), IR(4,:,2)/max(IR(4,:,2)), 'rs:')
title('40 degrees');
legend({'WT', 'AV'});
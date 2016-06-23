%% Calculate and Plot Bifurcations for cortical HH model
%==========================================================================
% This routine will recursively simulate steady state responses of a HH
% model adapted to resemble cortical neuronal function, for wild type
% sodium channels, and parameters derived from an epilepsy causing SCN1A
% mutation. Simulations will be performed at a range of input currents
% (defined as 'forward' variable below), with steady state values of one
% parameter step acting as initial conditions for the next parameter step.
% This is performed in two directions (forward and backward) in order to
% plot a bifurcation diagram of action potential generation. 
%
% The routine produces a figure that was the basis for Figure 6 in the
% publication below:
% 
% Peters C, Rosch RE, Hughes E, Ruben P (2016) Temperature-dependent 
% changes in neuronal dynamics in a patient with an SCN1A mutation and 
% hyperthermia induced seizures (under review)


% Manual definitions
%--------------------------------------------------------------------------
direction       = [1 -1];           % 1 = forward, -1 = backward
forward         = 0:1:90;           % defines range of I_stim parameter
n_exp           = {'WT37','WT40','AV37','AV40'};
exp_parameters  = [1 2 3 4];        % 1 = WT37 (norm), 2 = WT40, 3 = AV37, 3 = AV40;
count_e         = 1;

for e = exp_parameters

    
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
    
for d = direction
    
% Use variables
%--------------------------------------------------------------------------
i   = 0;                % counts number of model runs
cyc = length(exp_parameters);

if d == 1 
    steps = forward;
    disp('Forward');
else
    steps = flip(forward);
    disp('Back');
end

% Run model to steady state from random initial
%--------------------------------------------------------------------------
options     = odeset('InitialStep',0.0025,'MaxStep',0.05);
t_range     = [0 100];
x_ini       = [0 0 0 0];
[t,x]       = ode45(@(t,x)cort_variable_hh(t,x,params),t_range,x_ini,options);

% Run model iteratively across parameterspace
%--------------------------------------------------------------------------
for para = steps; 
    
    disp(length(steps)-i)
    params(8) = para;
    i = i+1;

    bif(1,i)    = para;
    options     = odeset('InitialStep',0.0025,'MaxStep',0.05);
    t_range     = [0 100];
    x_ini       = x(end,:);
    [t,x]       = ode45(@(t,x)cort_variable_hh(t,x,params),t_range,x_ini,options); 
    
    bif(2,i) = max(x(floor(3/5*end):end,1));
    bif(3,i) = min(x(floor(3/5*end):end,1));
    
    if bif(2,i)-bif(3,i)>0.5   
        [pks,locs] = findpeaks(x(floor(4/5*end):end));
        bif(4,i) = 1/mean(diff(t(locs)));
    else
        bif(4,i) = 0;
    end
end

% Draw figures: Top panel - bifurcation diagram, bottom panel - frequency
%==========================================================================
figure(2);
set(gcf, 'Position', [100 200 1200 300]);
hold on
subplot(3,cyc,[count_e cyc+count_e])

% Bifurcation diagram
%--------------------------------------------------------------------------
for i = 1:length(bif(1,:));
    if bif(2,i)-bif(3,i)<0.5    
        if d == 1
            scatter(bif(1,i), bif(2,i),'k.'); hold on
        else
            scatter(bif(1,i), bif(2,i),'r.'); hold on
        end
    else 
        if d == 1 
            scatter(bif(1,i), bif(2,i), 10, 'ko'); hold on
            scatter(bif(1,i), bif(3,i), 10, 'kx'); 
        else
            scatter(bif(1,i), bif(2,i), 10, 'ro'); hold on
            scatter(bif(1,i), bif(3,i), 10, 'rx');
        end
    end
end

ylabel('Membrane potential (mV)');
switch e 
    case 1, title('WT at 37 degrees');
    case 2, title('WT at 40 degrees');
    case 3, title('AV at 37 degrees');
    case 4, title('AV at 40 degrees');
end

% Frequency 
%--------------------------------------------------------------------------
subplot(3,cyc,2*cyc+count_e)
if d == 1                            
    plot(bif(1,:),bif(4,:), 'k');
    xlabel('Input Current');
    ylabel('Frequency');
else 
    plot(bif(1,:),bif(4,:), 'r');
    xlabel('Input Current');
    ylabel('Frequency');
end

end     % from for loop over directions
count_e = count_e+1;
clear x t x_ini bif
end     % from loop over experimental conditions
%% Calculate and Plot Bifurcations for cortical HH model
%==========================================================================
% This routine will simulate a Hodgkin-Huxley model adapted to capture
% dynamics of cortical neurons for different parameterisations: the
% wildtype sodium channels at 37, and 40 degrees; an SCN1A mutation causing
% epilepsy at 37, and 40 degrees; and all of these at two different input
% current strengths. 
%
% The routine produces a figure that was the basis for Figure 6A in the
% publication below:
% 
% Peters C, Rosch RE, Hughes E, Ruben P (2016) Temperature-dependent 
% changes in neuronal dynamics in a patient with an SCN1A mutation and 
% hyperthermia induced seizures (under review)

% Housekeeping
%--------------------------------------------------------------------------
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
exp_parameters  = [1 2 3 4];  % 1 = WT37 (norm), 2 = WT40, 3 = AV37, 3 = AV40;
plot_phase      = 0;
I_stim          = [0.2 45];         % I_stim    HH value = 0.200 nA
direction       = [1 -1];           % 1 = forward, -1 = backward
forward         = 0:2:90;           % defines range of I_stim parameter

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


for i = 1:2; 
% Run Model
%==========================================================================
params(8)   = I_stim(i);
options     = odeset('InitialStep',0.005,'MaxStep',0.05);
t_range     = [0 50];
x_ini       = [0 0 0 0];
[t,x]       = ode45(@(t,x)cort_variable_hh(t,x,params),t_range,x_ini,options);

% Plot time courses
% --------------------------------------------------------------------------
Figure1     = figure(1); 
xplot = linspace(1,length(t),5000);
xplot = floor(xplot);

if e <= 2,  sp_no = 0;
else        sp_no = 2; 
end

if rem(e,2) == 0,   col = 'r'; 
else                col = 'k';
end

subplot(1,4,sp_no+i);  
    plot(t(xplot),x(xplot,1), col); hold on
	title(n_exp{e});
    axis([0,20,-100,100]);
    set(gca, 'XTickLabel', [], 'Box', 'off');
    xlabel('time');
    if count > 1
        set(gca, 'YTickLabel', []);
    else
        ylabel('Membrane voltage in mV');  
    end
    
set(Figure1, 'Position', [100 100 1000 350]);  
subplot(1,4,1), title('WT low stimulation');
subplot(1,4,2), title('WT high stimulation');
subplot(1,4,3), title('AV low stimulation');
subplot(1,4,4), title('WT high stimulation'); legend({'37', '40'});
end

end
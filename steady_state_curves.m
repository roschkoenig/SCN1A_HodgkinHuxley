%% Simulate data from cortical parameterisation of HH
%==========================================================================
% This routine will simulate datapoints for the cortical parameterisation
% of the HH model according to Traub (1991) and fit the parameters of the
% steady state / voltage clamp formulation to these data in order to find
% the baseline for the further modelling steps
% 
% This script is the basis for Figure 5a and c in the publication below:
% 
% Peters C, Rosch RE, Hughes E, Ruben P (2016) Temperature-dependent 
% changes in neuronal dynamics in a patient with an SCN1A mutation and 
% hyperthermia induced seizures (under review)

%
% Setting up the HH equations
%--------------------------------------------------------------------------
clear all
syms a_m b_m a_h b_h V m_inf(V) h_inf(V) t_m(V) t_h(V)
V_t = -65;

a_m = -0.32 * ( V - V_t - 13 ) / (exp(-(V - V_t - 13) / 4) - 1);
a_h = 0.128 * exp( -( V - V_t - 17 )/18 );
b_m = 0.28 * ( V - V_t - 40 ) / ( exp( (V - V_t - 40)/5 )-1 );
b_h = 4 / ( 1 + exp( -(V-V_t-40)/5 ) );

m_inf(V) = a_m / (a_m + b_m);
h_inf(V) = a_h / (a_h + b_h);

t_m(V) = 1/(a_m + b_m);
t_h(V) = 1/(a_h + b_h);

% Calculating simulated datapoints
%--------------------------------------------------------------------------
i = 0;
for V = -80:3:20
    i = i+1;
    hh_cort(1,i) = V;
    hh_cort(2,i) = m_inf(V);
    hh_cort(3,i) = t_m(V);
    hh_cort(4,i) = h_inf(V);
    hh_cort(5,i) = t_h(V);
end

% Fitting steady state formulation to simulated datapoints
%--------------------------------------------------------------------------
clear V

V = hh_cort(1,:);
m = hh_cort(2,:);
h = hh_cort(4,:);

F_m   = @(x_m,xdata_m) 1 ./ (1 + exp(-(xdata_m-x_m(1))./x_m(2)));
F_h   = @(x_h,xdata_h) 1 ./ (1 + exp((xdata_h-x_h(1))./x_h(2)));
x0  = [-12,3];
[x_m, resnorm, ~, exitflag, output] = lsqcurvefit(F_m,x0,V,m);
[x_h, resnorm, ~, exitflag, output] = lsqcurvefit(F_h,x0,V,h);

% Plotting the estimated and simulated steady state curves
%--------------------------------------------------------------------------
for V = -80:3:20
    i = i+1;
    estimated(1,i) = V;
    estimated(2,i) = F_m(x_m,V);
    estimated(4,i) = F_h(x_h,V);
end

figure(1)
subplot(3,1,1)
    plot(estimated(1,:),estimated(2,:),'ro'); hold on
    plot(hh_cort(1,:), hh_cort(2,:), 'r');
    plot(estimated(1,:),estimated(4,:),'bo'); hold on
    plot(hh_cort(1,:), hh_cort(4,:), 'b');


% Plotting steady state curves for different experimental values
%==========================================================================
% This routine will take the experimentally derived voltage clamp values
% and plot the respective steady state curves for the different gating
% parameters of the Hodgkin Huxley Model


% Choosing parmeterisation
%--------------------------------------------------------------------------
n_exp = {'WT37','WT40','AV37','AV40'};
exp_parameters = 1:4; % 1 = WT37 (norm), 2 = WT40, 3 = AV37, 3 = AV40;


for e = exp_parameters
% Set according to literature
%--------------------------------------------------------------------------
params(1) = 0.010;      % C         HH
params(2) = 0.0000205; 	% g_L       Posposchil 2008, Fig 2
params(3) = 5;          % g_K       Traub 1991
params(4) = 56;         % g_Na      Posposchil 2008, Fig 2
params(5) = -90;        % E_K       Traub 1991
params(6) = -70.3;    	% E_L       Posposchil 2008, Fig 2
params(7) = 50;         % E_Na      Traub 1991
params(8) = .2;         % I_stim    HH
params(11) = -55;      	% V_t       arbitrary (guided by Posposchil)


% Set according to experimental values (with WT37 as baseline)
%--------------------------------------------------------------------------
switch e
    case 1  % standard HH formulation (WT37)
        temperature = 37 + 273;       
        dV2_m   = 0;
        dV2_h   = 0;
        t_off   = 0;
        m_z     =  3.6089;
        h_z     = -6.6210;
        
    case 2  % Experimental parameters for WT40
        temperature = 40 + 273;
        dV2_m   =  6.49;
        dV2_h   = -1.9;   
        t_off   = 1.9;
        m_z    	=  2.43*3.6089/2.14;    
        h_z   	= -3.79*6.6210/3.14;
        
    case 3  % Experimental parameters for AV37
        temperature = 37 + 273;
        dV2_m   =  5.94;
        dV2_h   =  2.4; 
        t_off   = -2.4;
        m_z    	=  2.58*3.6089/2.14;    
        h_z   	= -3.73*6.6210/3.14;
        
    case 4  % Experimental parameters for AV40
        temperature = 40 + 273;
        dV2_m   = 9.93;
        dV2_h   = 11.7;
        t_off   = -11.7;
        m_z    	= 2.1*3.6089/2.14;    
        h_z   	= -3.38*6.6210/3.14;      
end

s_m         = (0.0863 * temperature)/m_z;
s_h         = -(0.0863 * temperature)/h_z;

params(9) = x_m(1) + dV2_m;	% V_2m, baseline estimated above
params(10) = s_m;           % s_m

params(12) = x_h(1) + dV2_h; % V_2h, baseline estimated above
params(13) = -s_h;          % s_h
params(14) = t_off;   
%--------------------------------------------------------------------------

syms x1
m_inf = (1 / (1 + exp(-(x1-params(9))/s_m)));
h_inf = 1 / (1 + exp(-(x1-params(12))/-s_h));

figure(1);
if rem(e,2) == 0, col = 'r'; else col = 'k'; end
if e <= 2
    subplot(3,1,2)
    h = ezplot(h_inf, [-80 20]); hold on
    g = ezplot(m_inf, [-80 20]); hold on
    set(h, 'Color', col);
    set(g, 'Color', col);
    legend({'WT37 h','WT37 m','WT40 h','WT40 m'});
    title('WT channel');
else
    subplot(3,1,3)
    h = ezplot(h_inf, [-80 20]); hold on
    g = ezplot(m_inf, [-80 20]); hold on
    set(h, 'Color', col);
    set(g, 'Color', col);
    legend({'AV37 h','AV37 m','AV40 h','AV40 m'});
    title('AV mutation');
end

dis = figure(1);
set(dis, 'Position', [100 100 400 700]);
end
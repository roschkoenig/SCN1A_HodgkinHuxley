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
%         m_z  	= 2.14;
%         h_z  	= -3.14;
        
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

params(9) = -39.0 + dV2_m;	% V_2m
params(10) = s_m;           % s_m

params(12) = -43.3 + dV2_h; % V_2h
params(13) = -s_h;          % s_h

params(14) = t_off;         %


cols = jet(4);

syms x1
m_inf = 1 / (1 + exp(-(x1-params(9))/s_m));
h_inf = 1 / (1 + exp(-(x1-params(12))/-s_h));

figure(1);
h = ezplot(h_inf, [-80 20]); hold on
g = ezplot(m_inf, [-80 20]); hold on
set(h, 'Color', cols(e,:));
set(g, 'Color', cols(e,:));
legend({'WT37 h','WT37 m','WT40 h','WT40 m','AV37 h','AV37 m','AV40 h','AV40 m'});
end
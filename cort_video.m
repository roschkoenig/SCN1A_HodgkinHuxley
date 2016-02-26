% Choose experiments to be modelled
%--------------------------------------------------------------------------
n_exp           = {'WT37','WT40','AV37','AV40'};
e               = 1;  % 1 = WT37 (norm), 2 = WT40, 3 = AV37, 3 = AV40;
I_stim          = 0.2;         % I_stim    HH value = 0.200 nA

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


% Run Model
%==========================================================================
options     = odeset('InitialStep',0.0005,'MaxStep',0.005);
t_range     = [0 50];
x_ini       = [0 0 0 0];
[t,x]       = ode45(@(t,x)cort_variable_hh(t,x,params),t_range,x_ini,options);

%% Plot Parameters and Phase Space
% --------------------------------------------------------------------------
Figure1     = figure;
set(gcf,'color','w');
set(Figure1, 'PaperUnits', 'centimeters', 'PaperSize',[21 15],...
    'PaperPosition',[0,0,6.5,4.5]);   

xplot = linspace(1, length(t),1000);
xplot = floor(xplot);
subplot(2,3,1) 
    plot(t,x(:,1), 'LineWidth', 2); hold on;
    timetrace = gca;
    set(timetrace, 'XLim', [-5 0])
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', [], 'Box', 'off');

subplot(2,3,4)   
    plot(t,x(:,2:end), 'LineWidth', 2);
    parametertrace = gca;
    set(parametertrace, 'XLim', [-5 0])
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', [], 'Box', 'off');
    
pause
last = xplot(1); 
for i = xplot(2:1:length(xplot))
    
    set(timetrace, 'XLim', [-5 0]+t(i));
    set(parametertrace, 'XLim', [-5 0]+t(i));

    subplot(2,3,[2:3,5:6])
        delete(last_dot)
        last_dot = scatter3(x(i,2), x(i,3), x(i,4), 100, 'r', 'filled');
        scatter3(x(i,2), x(i,3), x(i,4), 'k', 'filled');
        axis([0 1 0 1 0 1]);
        view(-15.5,40)
        hold on;
    pause(0.00001);
end
    
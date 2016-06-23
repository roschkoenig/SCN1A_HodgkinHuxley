function dxdt = cort_variable_hh(t,x,params, p1, p2)
dxdt=zeros(4,1);

% Constant values
%--------------------------------------------------------------------------
C       = params(1); 
g_L     = params(2); 
g_K     = params(3); 
g_Na    = params(4); 
E_K     = params(5); 
E_L     = params(6); 
E_Na    = params(7); 

I_stim  = params(8);

V_2m    = params(9);
s_m     = params(10);
V_t     = params(11);

V_2h    = params(12);
s_h     = params(13);
t_off     = params(14);

% Pulse definition
%==========================================================================
if t < p1(1)                         % Initial condition
    vol = -130;
elseif t > p1(1) && t < p1(2)     % Pulse 1
    vol = 0;
elseif t > p2(1) && t < p2(2)     % Pulse 2 
    vol = 0;
else                                    % Recovery window
    vol = -90;
end

% Parameter equations
%==========================================================================
% Opening probabilities 
%--------------------------------------------------------------------------
aV_m = -0.32 * ( vol - V_t - 13 ) / (exp(-(vol - V_t - 13) / 4) - 1);
aV_h = 0.128 * exp( -( vol - V_t - 17 + t_off )/18 );

% Closing probabilities
%--------------------------------------------------------------------------
bV_m = 0.28 * ( vol - V_t - 40 ) / ( exp( (vol - V_t - 40)/5 )-1 );
bV_h = 4 / ( 1 + exp( -(vol-V_t-40+t_off)/5 ) );

% Steady state values
%--------------------------------------------------------------------------
m_inf = 1 / (1 + exp(-(vol-V_2m)/s_m));
h_inf = 1 / (1 + exp(-(vol-V_2h)/-s_h));
t_h   = 1/(aV_h + bV_h);
t_m   = 1/(aV_m + bV_m);

% Differential equations
%--------------------------------------------------------------------------
% dxdt(1) = Na channel conductance current
% dxdt(3) = m sodium channel activation gate
% dxdt(4) = h sodium channel inactivation gate

dxdt(1) = (m_inf - x(1)) / t_m;             % m
dxdt(2) = (h_inf - x(2)) / t_h;             % h

end

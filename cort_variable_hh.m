function dxdt = cort_variable_hh(t,x,params)
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

% Parameter equations
%==========================================================================
% Opening probabilities 
%--------------------------------------------------------------------------
aV_n = -0.032 * ( x(1) - V_t - 15 ) / ( exp( -(x(1)-V_t-15)/5)-1 );
aV_m = -0.32 * ( x(1) - V_t - 13 ) / (exp(-(x(1) - V_t - 13) / 4) - 1);
aV_h = 0.128 * exp( -( x(1) - V_t - 17 + t_off )/18 );

% Closing probabilities
%--------------------------------------------------------------------------
bV_n = 0.5 * exp( -(x(1) - V_t - 10)/40);
bV_m = 0.28 * ( x(1) - V_t - 40 ) / ( exp( (x(1) - V_t - 40)/5 )-1 );
bV_h = 4 / ( 1 + exp( -(x(1)-V_t-40+t_off)/5 ) );

% Steady state values
%--------------------------------------------------------------------------
m_inf = 1 / (1 + exp(-(x(1)-V_2m)/s_m));
h_inf = 1 / (1 + exp(-(x(1)-V_2h)/-s_h));
t_h   = 1/(aV_h + bV_h);
t_m   = 1/(aV_m + bV_m);

% Differential equations
%--------------------------------------------------------------------------
% dxdt(1) = V
% dxdt(2) = n potassium channel activation gate
% dxdt(3) = m sodium channel activation gate
% dxdt(4) = h sodium channel inactivation gate

dxdt(1) = ( -g_Na * x(3)^3 * x(4) * (x(1) - E_Na)  ...
            -g_K * x(2)^4 * (x(1) - E_K)  ...
            -g_L * (x(1) - E_L) + I_stim ) / C;  
        
dxdt(2) = aV_n * (1-x(2)) - bV_n * x(2);    % n
dxdt(3) = (m_inf - x(3)) / t_m;             % m
dxdt(4) = (h_inf - x(4)) / t_h;             % h

end

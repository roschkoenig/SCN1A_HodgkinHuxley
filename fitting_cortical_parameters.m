% Simulate data from cortical parameterisation of HH
%==========================================================================
% This routine will simulate datapoints for the cortical parameterisation
% of the HH model according to Traub (1991) and fit the parameters of the
% steady state / voltage clamp formulation to these data in order to find
% the baseline for the further modelling steps

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

plot(estimated(1,:),estimated(2,:),'ro'); hold on
plot(hh_cort(1,:), hh_cort(2,:), 'r');
plot(estimated(1,:),estimated(4,:),'bo'); hold on
plot(hh_cort(1,:), hh_cort(4,:), 'b');
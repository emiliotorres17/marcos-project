close all
clear all
clc
%-------------------------------------------------------------------------%
% Defining global variables                                               %
%-------------------------------------------------------------------------%
global k 
global mu_mix 
global mu_deep 
global lambda 
global n 
global gamma
%-------------------------------------------------------------------------%
% Defining domain variables                                               %
%-------------------------------------------------------------------------%
k       = 5.35*3.154e+07;
mu_deep = 6.307e+09;
mu_mix  = 3.154e+08;
n       = k*log(2);
gamma   = 1.2*3.154e+07;
lambda  = 1.2*3.154e+07;
%-------------------------------------------------------------------------%
% solving ODE                                                             %
%-------------------------------------------------------------------------%
[t, Ta] = ode45('func', [0 1000], [0 0]);
T_mix   = Ta(:,1);
T_deep  = Ta(:,2);
%-------------------------------------------------------------------------%
% Plotting                                                                %
%-------------------------------------------------------------------------%
hold on
plot(t, T_mix, 'r', 'linewidth', 1.5)
plot(t, T_deep, 'b', 'linewidth', 1.5)
legend('$\Delta T_{mix}$', '$\Delta T_{deep}$', 'interpreter', 'latex')
xlabel('Time [yrs]', 'interpreter', 'latex')
ylabel('Temperature', 'interpreter', 'latex')
grid('on')

function dT = my_ode(t,T)
    %mu_deep = 6.307e+09;
    %mu_mix  = 3.154e+08;
    %lambda  = 1.2*3.154e+07;
    %gamma   = 1.2*3.154e+07;
    %k       = 5.35*3.154e+07;
    %n       = k*log(2);
    dT_mix      = (1.0/mu_mix)*(-lambda*T(1) - ...
                    gamma*(T(1) -T(2)) + n);
    dT_deep     = (1.0/mu_deep)*(gamma*(T(1) - T(2)));

    dT           = [dT_mix;  dT_deep];
end 

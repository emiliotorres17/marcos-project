close all
clear all
clc
%-------------------------------------------------------------------------%
% Defining domain variables                                               %
%-------------------------------------------------------------------------%
k       = 5.35*3.154e+07;
lambda  = 1.2*3.154e+07;
mu_deep = 6.307e+09;
mu_mix  = 3.154e+08;
gamma   = 1.2*3.154e+07;
n       = k*log(2);
%-------------------------------------------------------------------------%
% Defining time variables                                                 %
%-------------------------------------------------------------------------%
dt      = 0.5;
t       = 0:dt:1000;
T_mix   = zeros(1,length(t)); 
T_deep  = zeros(1,length(t)); 
%-------------------------------------------------------------------------%
% Time loop                                                               %
%-------------------------------------------------------------------------%
for i = 1:length(t)-1        
    dT_mix      = (1.0/mu_mix)*(-lambda*T_mix(i) - ...
                    gamma*(T_mix(i) -T_deep(i)) + n);
    dT_deep     = (1.0/mu_deep)*(gamma*(T_mix(i) - T_deep(i)));
    T_mix(i+1)  = T_mix(i) + dt*dT_mix;
    T_deep(i+1) = T_deep(i) + dt*dT_deep;
end
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

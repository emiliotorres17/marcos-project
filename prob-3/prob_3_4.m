clear all
close all
clc
%-------------------------------------------------------------------------%
% Defining global variables                                               %
%-------------------------------------------------------------------------%
global k
global lambda
global mu_deep
global mu_mix
global gamma
global ka
global kd
global del_d
global AM
global OM
global kH
global k1
global k2
global alk
global del_a
%global E
%-------------------------------------------------------------------------%
% Reading in data                                                         %
%-------------------------------------------------------------------------%
Em                  = xlsread('../data/EmissionsProj.xlsx', 'B2:B265');
El                  = xlsread('../data/EmissionsProj.xlsx', 'C2:C265');
%Em                  = (Em+El)*10^15/12.0;
Em                  = (Em)*10^15/12.0;
emt                 = linspace(1751, 2014, length(Em));
land_emissions      = xlsread('../data/EmissionsProj.xlsx', 'C2:C265');
emissions_yrs       = xlsread('../data/EmissionsProj.xlsx', 'A2:A265');
temp_data           = xlsread('../data/TempProj.xlsx', 'B2:B170');
temp_yrs            = xlsread('../data/TempProj.xlsx', 'A2:A170');
CO2_data            = xlsread('../data/CO2_proj.xlsx', 'B2:B106');
CO2_yrs             = xlsread('../data/CO2_proj.xlsx', 'A2:A106');
%-------------------------------------------------------------------------%
% Temperature variables                                                   %
%-------------------------------------------------------------------------%
k       = 5.35*3.154e+07;
%lambda  = 1.2*3.154e+07;
mu_deep = 6.307e+09;
mu_mix  = 3.154e+08;
gamma   = 1.2*3.154e+07;
%-------------------------------------------------------------------------%
% CO2 variables                                                           %
%-------------------------------------------------------------------------%
ka      = 0.2;
kd      = 0.05;
del_d   = 50.0;
AM      = 1.77e+20;
OM      = 7.8e+22;
kH      = 1.23e+3;
k1      = 1.44e-08;
k2      = 8.154e-12;
alk     = 767e+15/12.0;
del_a   = OM/(AM*(1+del_d));

Qa      = 590*10^15/12;     
Qu      = 713*10^15/12;     
Ql      = 35658*10^15/12.0; 
%-------------------------------------------------------------------------%
% solving for CO2                                                         %
%-------------------------------------------------------------------------%
Qint    = [Qa; Qu; Ql]; 
options        = odeset('MaxStep', 1);
[t, out]       = ode45(@(t,y) func2(t, y, emt, Em), [1751 2014], Qint);
hold on
%plot(t, temp)
Qa  = out(:,1)./AM.*10^6;
plot(t, Qa, 'r', 'linewidth', 1.5)
plot(CO2_yrs, CO2_data, 'b', 'linewidth', 1.5)
legend('Model prediction', 'Historical CO$_{2}$', 'interpreter', 'latex', ...
        'location' , 'Northwest')
xlabel('Time [yrs]', 'interpreter', 'latex')
ylabel('Atmospheric CO$_{2}$', 'interpreter', 'latex' )
grid('on')
xlim([1751, 2014])

%-------------------------------------------------------------------------%
% solving for temperature                                                 %
%-------------------------------------------------------------------------%
C           = Qa;
C0          = 277.0;
n           = k*log(C/C0); 
lambda      = gamma;
nt          = linspace(1751, 2014, length(n));
Tint        = [0; 0];
[t, out]    = ode45(@(t,y) func3(t, y, nt, n), [1751 2014], Tint);
Tmix        = out(:,1);
figure          
hold on
plot(t, Tmix, 'r', 'linewidth', 1.5)
plot(temp_yrs, temp_data, 'b', 'linewidth', 1.5)
legend('Model Prediction', 'Historical Temperature Anomaly', ...
            'interpreter', 'latex', 'location' , 'Northwest')
xlabel('Time [yrs]', 'interpreter', 'latex')
ylabel('Temperature Anomaly [K]', 'interpreter', 'latex' )
grid('on')
xlim([1751, 2014])

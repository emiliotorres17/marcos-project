close all
clear all
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
%-------------------------------------------------------------------------%
% Plotting flags                                                          %
%-------------------------------------------------------------------------%
E_flag      = 1;
CO2_flag    = 1;
T_flag      = 1;
n_flag      = 1;
%-------------------------------------------------------------------------%
% Defining emissions                                                      %
%-------------------------------------------------------------------------%
emt     = 2020:1:2300;
Em      = zeros(length(emt), 1);
percent = 100.0
E_ramp  = 1;
if E_ramp == 1
    for i =1:length(Em)
        if emt(i) <= 2050
            Em(i) = -2/5*emt(i) + 820;
        else
            Em(i) = 0.0;
        end
        if Em(i) < 0.0
            Em(i) = 0.0;
        end
    end
else
    Em = 12.0 + Em;
end
%-------------------------------------------------------------------------%
% Plotting E(t)                                                           %
%-------------------------------------------------------------------------%
h(1) = subplot(3,2,1);
plot(emt, Em, 'r', 'linewidth', 1.5);
legend('Emissions', 'interpreter', 'latex', 'location' , 'Northeast')
xlabel('Time [yrs]', 'interpreter', 'latex')
ylabel('Emissions [GtC]', 'interpreter', 'latex' )
grid('on')
xlim([2020, 2300])
Em  = (Em)*10^15/12.0;
%-------------------------------------------------------------------------%
% Reading in data                                                         %
%-------------------------------------------------------------------------%
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
Qa      = 870*10^15/12;     
Qu      = 730*10^15/12;     
Ql      = 35707*10^15/12.0; 
%-------------------------------------------------------------------------%
% solving for CO2                                                         %
%-------------------------------------------------------------------------%
Qint            = [Qa; Qu; Ql]; 
options         = odeset('MaxStep', 1);
[t_co2, out]    = ode45(@(t,y) func2(t, y, emt, Em), [2020 2300], Qint, ...
                        options);
%-------------------------------------------------------------------------%
% Plotting CO2                                                            %
%-------------------------------------------------------------------------%
Qa  = out(:,1)./AM.*10^6;
h(2) = subplot(3,2,2);
plot(t_co2, Qa, 'r', 'linewidth', 1.5);
legend('Model prediction', 'interpreter', 'latex', ...
        'location' , 'best')
xlabel('Time [yrs]', 'interpreter', 'latex')
ylabel('Atmospheric CO$_{2}$ [ppm]', 'interpreter', 'latex' )
grid('on')
xlim([2020, 2300])
%-------------------------------------------------------------------------%
% Setting n(t)                                                            %
%-------------------------------------------------------------------------%
n_ramp      = 1;
no_albedo   = 1;
C           = Qa;
C0          = 277.0;
n0          = k*log(C(1)/C0);
nt          = linspace(2020, 2300, length(C))';
n           = zeros(length(nt), 1);
if n_ramp == 1
    for i = 1:length(n)
        n(i) = -n0/100.0*nt(i) + 2050*n0/100.0;
        if n(i) < 0.0
            n(i) = 0.0;
        end
    end
else
    n = 0.0*n;
end
if no_albedo == 1
    n   = k*log(C/C0);
    disp('test');
end
nstore  = n;
%-------------------------------------------------------------------------%
% solving for temperature                                                 %
%-------------------------------------------------------------------------%
lambda      = gamma;
Tint        = [0.94; 0.25];
options     = odeset('MaxStep', 1);
[t_T, out_T]= ode45(@(t,y) func3(t, y, nt, n), [2020 2300], Tint, options);
Tmix        = out_T(:,1);
disp(max(Tmix));
%-------------------------------------------------------------------------%
% Plotting temperature                                                    %
%-------------------------------------------------------------------------%
h(3) = subplot(3,2,3);
hold on
plot(t_T, Tmix, 'r', 'linewidth', 1.5);
h1  = plot([2020, 2300], [1.5, 1.5], 'b--', 'linewidth', 1.5);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Model Prediction',  'interpreter', 'latex', 'location' , 'best');
xlabel('Time [yrs]', 'interpreter', 'latex');
ylabel('Temperature Anomaly [K]', 'interpreter', 'latex' );
grid('on')
xlim([2020, 2300])
%-------------------------------------------------------------------------%
% Plotting N(t)                                                           %
%-------------------------------------------------------------------------%
h(4) = subplot(3,2,4);
plot(t_co2, nstore, 'r', 'linewidth', 1.5);
legend('Radiative Force',  'interpreter', 'latex', 'location' , 'Southeast');
xlabel('Time [yrs]', 'interpreter', 'latex');
ylabel('Radiative Force [Wm$^{-2}$]', 'interpreter', 'latex' );
grid('on')
xlim([2020, 2300])
%---------------------------------------------------------------------%
% Calculating H                                                       %
%---------------------------------------------------------------------%
Qu  = out(:,2); 
ph  = zeros(length(t_co2),1);
for i = 1:length(Qu) 
    qu  = Qu(i);
    b   = k1*(1-qu/alk);
    c   = k1*k2*(1 - (2.0*qu)/alk);
    a   = 1.0;
    H   = (-b - (b^2 - 4*a*c)^(0.5))/(2.0*a);
    if H < 0.0
        H = (-b + (b^2 - 4*a*c)^(0.5))/(2.0*a);
    end
    ph(i)   = - log10(H*(1000/18));
end
h(5) = subplot(3,2,5);
plot(t_co2, ph, 'r', 'linewidth', 1.5);
pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)])
legend('Ocean pH',  'interpreter', 'latex', 'location' , 'Southeast');
xlabel('Time [yrs]', 'interpreter', 'latex');
ylabel('Ocean pH', 'interpreter', 'latex' );
grid('on')
xlim([2020, 2300])

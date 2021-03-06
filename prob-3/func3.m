function [dT] = func3(t,T, nt, n)
%-------------------------------------------------------------------------%
% Defining global variables                                               %
%-------------------------------------------------------------------------%
    global k 
    global mu_mix 
    global mu_deep 
    global lambda 
    global gamma
    n = interp1(nt, n, t);
    dT_mix      = (1.0/mu_mix)*(-lambda*T(1) - ...
                    gamma*(T(1) -T(2)) + n);
    dT_deep     = (1.0/mu_deep)*(gamma*(T(1) - T(2)));

    dT           = [dT_mix;  dT_deep];
end 

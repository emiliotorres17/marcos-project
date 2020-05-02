function dQ = func2(t,Q, emt, Em)
%-------------------------------------------------------------------------%
% Calling global variables                                                %
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
global i
%---------------------------------------------------------------------%
% Calculating H                                                       %
%---------------------------------------------------------------------%
qu  = Q(2);
b   = k1*(1-qu/alk);
c   = k1*k2*(1 - (2.0*qu)/alk);
a   = 1.0;
H   = (-b - (b^2 - 4*a*c)^(0.5))/(2.0*a);
if H < 0.0
    H = (-b + (b^2 - 4*a*c)^(0.5))/(2.0*a);
end
%---------------------------------------------------------------------%
% Calculating lambda                                                  %
%---------------------------------------------------------------------%
lambda  = 1 + k1/H + (k1*k2)/H^2;

Em = interp1(emt, Em, t);
%---------------------------------------------------------------------%
% Derivatives                                                         %
%---------------------------------------------------------------------%
d_Qa    = -ka*Q(1) + ka*kH/(del_a*lambda)*Q(2) + Em;
d_Qu    = ka*Q(1) - ka*kH/(del_a*lambda)*Q(2) - kd*Q(2) + kd/del_d *Q(3);
d_Ql    = kd*Q(2) -  kd/del_d*Q(3);

dQ      = [d_Qa; d_Qu; d_Ql]; 


end

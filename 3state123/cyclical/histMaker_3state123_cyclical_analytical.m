function [Peq] = histMaker_3state123_cyclical_analytical(t12,t13,t21,t23,t31,A1,A2,A3)
% function [Peq] = histMaker_3state123_cyclical_analytical(k12,k21,k23,k32,k31)
% k13 = k31*k23*k12/(k32*k21); % This relation forces detailed balance to be satisfied.
%--------------------------------------------------------------------------
% Set the rates
%--------------------------------------------------------------------------
k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;
% Detailed balance condition:
k32 = k12*k23*k31/(k13*k21);

% Define equilibrium values of populations
p1_eq = 1/(1 + (k12 + k13)/k21 + (1 - k31/k21)*((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21));
p3_eq = ((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21)*p1_eq;
p2_eq = ((k12 + k13)*p1_eq - k31*p3_eq)/k21;

Peq = [p1_eq; p2_eq; p3_eq];
end


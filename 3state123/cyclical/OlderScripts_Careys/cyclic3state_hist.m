function [p1_eq,p2_eq,p3_eq] = cyclic3state_hist(A1,A2,A3,k12,k21,k23,k32,k31)
%This gives the same output as the new numerical methods

k13 = k31*k23*k12/(k32*k21); % This relation forces detailed balance to be satisfied.

% Define equilibrium values of populations
p1_eq = 1/(1 + (k12 + k13)/k21 + (1 - k31/k21)*((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21));
p3_eq = ((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21)*p1_eq;
p2_eq = ((k12 + k13)*p1_eq - k31*p3_eq)/k21;

end
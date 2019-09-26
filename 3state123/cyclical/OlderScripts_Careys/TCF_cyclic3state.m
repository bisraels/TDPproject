function tcf = TCF_cyclic3state(time,A0,A1,A2,k01,k10,k12,k21,k20)
% A = [A0,A1,A2];
% A0 = A0 - mean(A);
% A1 = A1 - mean(A);
% A2 = A2 - mean(A);

k02 = k20*k12*k01/(k21*k10); % This relation forces detailed balance to be satisfied.

% Define eigenvalues lam1 and lam2
c1 = 0.5*(k01+k10+k12+k21+k02+k20);
c2 = 0.5*sqrt(k01^2 + 2*k01*k02 + 2*k01*k10 - 2*k01*k12 - 2*k01*k20 - 2*k01*k21 + k02^2 - 2*k02*k10 - 2*k02*k12 + 2*k02*k20 - k02*k21 + k10^2 + 2*k10*k12 - 2*k10*k20 - 2*k10*k21 + k12^2 - 2*k12*k20 + 2*k12*k21 + k20^2 + 2*k20*k21 + k21^2);
lam1 = c1 + c2;
lam2 = c1 - c2;

% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
a = (k12 + k20 + k21 - c1 - c2)/(k02 - k12);
b = (c1 + c2 - k02 - k20 - k21)/(k02 - k12);
c = 1;

% Define eigenvector components corresponding to lam2: v2 = [x,y,z].
x = (k12 + k20 + k21 - c1 + c2)/(k02 - k12);
y = (c1 - c2 - k02 - k20 - k21)/(k02 - k12);
z = 1;

% Define equilibrium values of populations
p0_eq = 1/(1 + (k01 + k02)/k10 + (1 - k20/k10)*((k10 + k12)*(k01 + k02) - k01*k10)/((k10 + k12)*k20 + k21*k10));
p2_eq = ((k10 + k12)*(k01 + k02) - k01*k10)/((k10 + k12)*k20 + k21*k10)*p0_eq;
p1_eq = ((k01 + k02)*p0_eq - k20*p2_eq)/k10;

% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_2 = (-p1_eq + b/a*p0_eq)/(y - b/a*x);
n_1 = (1 - p1_eq + b/a*p0_eq)/(y - b/a*x);
n_0 = (-p1_eq - b/a*(1 - p0_eq))/(y - b/a*x);
m_2 = (-p0_eq - n_2*x)/a;
m_1 = (-p0_eq - n_1*x)/a;
m_0 = (1 - p0_eq - n_0*x)/a;

% Define populations as a function of time for different initial values.
% p1_0 population 1 as a function of time, with the population initially in
% state 0.
p2_2 = p2_eq + m_2*c*exp(-lam1*time) + n_2*z*exp(-lam2*time);
p2_1 = p2_eq + m_1*c*exp(-lam1*time) + n_1*z*exp(-lam2*time);
p2_0 = p2_eq + m_0*c*exp(-lam1*time) + n_0*z*exp(-lam2*time);
p1_2 = p1_eq + m_2*b*exp(-lam1*time) + n_2*y*exp(-lam2*time);
p1_1 = p1_eq + m_1*b*exp(-lam1*time) + n_1*y*exp(-lam2*time);
p1_0 = p1_eq + m_0*b*exp(-lam1*time) + n_0*y*exp(-lam2*time);
p0_2 = p0_eq + m_2*a*exp(-lam1*time) + n_2*x*exp(-lam2*time);
p0_1 = p0_eq + m_1*a*exp(-lam1*time) + n_1*x*exp(-lam2*time);
p0_0 = p0_eq + m_0*a*exp(-lam1*time) + n_0*x*exp(-lam2*time);

% Subtract mean values
Amean = p0_eq*A0 + p1_eq*A1 + p2_eq*A2;
A0 = A0 - Amean;
A1 = A1 - Amean;
A2 = A2 - Amean;

% Calculate TCF
tcf = p2_eq*A2*(A2*p2_2 + A1*p1_2 + A0*p0_2) +...
    p1_eq*A1*(A2*p2_1 + A1*p1_1 + A0*p0_1) +...
    p0_eq*A0*(A2*p2_0 + A1*p1_0 + A0*p0_0);
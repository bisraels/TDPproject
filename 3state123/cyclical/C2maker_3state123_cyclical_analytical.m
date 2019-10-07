% function tcf = TCF_cyclic3state(time,A0,A1,A2,k01,k10,k12,k21,k20)
% function tcf = C2maker_3state123_cyclical_analytical(time,A1,A2,A3,k12,k21,k23,k32,k31)
function tcf = C2maker_3state123_cyclical_analytical(t12,t13,t21,t23,t31,A1,A2,A3,time)
% Now is a mean subtracted 2-point TCF. Works the same as TCF_cyclic3state

% note: this messes it up:
%t12 = 0.000825, t13 = 0.006767, t21 = 0.000615, t23 = 0.006725, t31 = 0.003205, t32 = 0.004270  A1 = 0.737326, A2 = 0.518160, A3 = 0.315085,
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

% Define eigenvalues lam1 and lam2
c1 = 0.5*(k12+k21+k23+k32+k13+k31);
% c2 = 0.5*sqrt(k12^2 + 2*k12*k13 + 2*k12*k21 - 2*k12*k23 - 2*k12*k31 - 2*k12*k32 + k13^2 - 2*k13*k21 - 2*k13*k23 + 2*k13*k31 - k13*k32 + k21^2 + 2*k21*k23 - 2*k21*k31 - 2*k21*k32 + k23^2 - 2*k23*k31 + 2*k23*k32 + k31^2 + 2*k31*k32 + k32^2);
c2 = 0.5*sqrt(k12^2 + 2*k12*k13 + 2*k12*k21 - 2*k12*k23 - 2*k12*k31 - 2*k12*k32 + k13^2 - 2*k13*k21 - 2*k13*k23 + 2*k13*k31 - 2*k13*k32 + k21^2 + 2*k21*k23 - 2*k21*k31 - 2*k21*k32 + k23^2 - 2*k23*k31 + 2*k23*k32 + k31^2 + 2*k31*k32 + k32^2);

lam1 = c1 + c2;
lam2 = c1 - c2;

% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
a = (k23 + k31 + k32 - c1 - c2)/(k13 - k23);
b = (c1 + c2 - k13 - k31 - k32)/(k13 - k23);
c = 1;

% Define eigenvector components corresponding to lam2: v2 = [x,y,z].
x = (k23 + k31 + k32 - c1 + c2)/(k13 - k23);
y = (c1 - c2 - k13 - k31 - k32)/(k13 - k23);
z = 1;

% Define equilibrium values of populations
p1_eq = 1/(1 + (k12 + k13)/k21 + (1 - k31/k21)*((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21));
p3_eq = ((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21)*p1_eq;
p2_eq = ((k12 + k13)*p1_eq - k31*p3_eq)/k21;

% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_3 = (-p2_eq + b/a*p1_eq)/(y - b/a*x);
n_2 = (1 - p2_eq + b/a*p1_eq)/(y - b/a*x);
n_1 = (-p2_eq - b/a*(1 - p1_eq))/(y - b/a*x);
m_3 = (-p1_eq - n_3*x)/a;
m_2 = (-p1_eq - n_2*x)/a;
m_1 = (1 - p1_eq - n_1*x)/a;

% Define populations as a function of time for different initial values.
% p1_0 population 1 as a function of time, with the population initially in
% state 0.
p3_3 = p3_eq + m_3*c*exp(-lam1*time) + n_3*z*exp(-lam2*time);
p3_2 = p3_eq + m_2*c*exp(-lam1*time) + n_2*z*exp(-lam2*time);
p3_1 = p3_eq + m_1*c*exp(-lam1*time) + n_1*z*exp(-lam2*time);
p2_3 = p2_eq + m_3*b*exp(-lam1*time) + n_3*y*exp(-lam2*time);
p2_2 = p2_eq + m_2*b*exp(-lam1*time) + n_2*y*exp(-lam2*time);
p2_1 = p2_eq + m_1*b*exp(-lam1*time) + n_1*y*exp(-lam2*time);
p1_3 = p1_eq + m_3*a*exp(-lam1*time) + n_3*x*exp(-lam2*time);
p1_2 = p1_eq + m_2*a*exp(-lam1*time) + n_2*x*exp(-lam2*time);
p1_1 = p1_eq + m_1*a*exp(-lam1*time) + n_1*x*exp(-lam2*time);

% Subtract mean values
Amean = p1_eq*A1 + p2_eq*A2 + p3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

% Calculate TCF
tcf = p3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
    p2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
    p1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
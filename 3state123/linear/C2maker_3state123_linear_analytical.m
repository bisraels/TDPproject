function C2_sim = C2maker_3state123_linear_analytical(t12,t21,t23,t32,A1,A2,A3,time)
%MODEL: 1 <--> 2 <--> 3
%--------------------------------------------------------------------------
% Set the rates
%--------------------------------------------------------------------------
k12 = 1/t12;
k21 = 1/t21;
k23 = 1/t23;
k32 = 1/t32;

% Define eigenvalues lam1 and lam2\  
c1 = 0.5*(k12+k21+k23+k32);
c2 = 0.5*sqrt(k12^2 + 2*k12*k21 - 2*k12*k23 - 2*k12*k32 + k21^2 + 2*k21*k23 - 2*k21*k32 + k23^2 + 2*k23*k32 + k32^2);
lam1 = c1 + c2;
lam2 = c1 - c2;

% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
a = (c1 + c2)/k21 - (k12 + k21)/k21;
b = k12/k21 - (c1 + c2)/k21;
c = 1;

% Define eigenvector components corresponding to lam2: v2 = [x,y,z].
x = (c1 - c2)/k21 - (k12 + k21)/k21;
y = k12/k21 - (c1 - c2)/k21;
z = 1;

% Define equilibrium values of populations
p2_eq = 1/(k21/k12 + 1 + k23/k32);
p1_eq = k21/k12*p2_eq;
p3_eq = k23/k32*p2_eq;

% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_3 = (-p2_eq - b/a*(1 - p3_eq))/(y - b/a*x);
n_2 = (1 - p2_eq + b/a*p3_eq)/(y - b/a*x);
n_1 = (-p2_eq + b/a*p3_eq)/(y - b/a*x);
m_3 = (1 - p3_eq - n_3*x)/a;
m_2 = (-p3_eq - n_2*x)/a;
m_1 = (-p3_eq - n_1*x)/a;

% Define populations as a function of time for different initial values.
% p2_1 population 2 as a function of time, with the population initially in
% state 1.
p3_3 = p3_eq + m_3*a*exp(-lam1*time) + n_3*x*exp(-lam2*time);
p3_2 = p3_eq + m_2*a*exp(-lam1*time) + n_2*x*exp(-lam2*time);
p3_1 = p3_eq + m_1*a*exp(-lam1*time) + n_1*x*exp(-lam2*time);
p2_3 = p2_eq + m_3*b*exp(-lam1*time) + n_3*y*exp(-lam2*time);
p2_2 = p2_eq + m_2*b*exp(-lam1*time) + n_2*y*exp(-lam2*time);
p2_1 = p2_eq + m_1*b*exp(-lam1*time) + n_1*y*exp(-lam2*time);
p1_3 = p1_eq + m_3*c*exp(-lam1*time) + n_3*z*exp(-lam2*time);
p1_2 = p1_eq + m_2*c*exp(-lam1*time) + n_2*z*exp(-lam2*time);
p1_3 = p1_eq + m_1*c*exp(-lam1*time) + n_1*z*exp(-lam2*time);

% Subtract mean values
Amean = p1_eq*A1 + p2_eq*A2 + p3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

% Calculate TCF
C2_sim = p3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) + p2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) + p1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_3);

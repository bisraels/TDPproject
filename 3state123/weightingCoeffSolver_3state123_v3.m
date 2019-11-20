% Define the eigenvector components
syms v1_1 v1_2 v1_3 % compoents of zeroth eigenvector
syms v2_1 v2_2 v2_3 %Components of the first eigenvector
syms v3_1 v3_2 v3_3 %Components of the second eigenvector

% Define expansion coefficients c's for our constants
syms c1_1 c1_2 c1_3  %Components of the zeroth coefficient
syms c2_1 c2_2 c2_3 %Components of the first coefficient
syms c3_1 c3_2 c3_3 %Components of the second coefficient

% Define the eqilibrium probabilities
syms p1_eq p2_eq p3_eq %Equilibrium Populations

% probability from i to j
% pji(t=0) = c1_i*v1_j + c2_i*v2_j + c3_i*v3_j
% P i --> 1   (Final state 1)
eqn1 = 1 == c1_1*v1_1 + c2_1*v2_1 + c3_1*v3_1;   % P11
eqn2 = 0 == c1_2*v1_1 + c2_2*v2_1 + c3_2*v3_1;   % P12
eqn3 = 0 == c1_3*v1_1 + c2_3*v2_1 + c3_3*v3_1;   % P13

% P i --> 2   (Final state 2)
eqn4 = 0 == c1_1*v1_2 + c2_1*v2_2 + c3_1*v3_2;   % P21
eqn5 = 1 == c1_2*v1_2 + c2_2*v2_2 + c3_2*v3_2;   % P22
eqn6 = 0 == c1_3*v1_2 + c2_3*v2_2 + c3_3*v3_2;   % P23

% P i --> 3   (Final state 3)                      
eqn7 = 0 == c1_1*v1_3 + c2_1*v2_3 + c3_1*v3_3;   % P31
eqn8 = 0 == c1_2*v1_3 + c2_2*v2_3 + c3_2*v3_3;   % P32
eqn9 = 1 == c1_3*v1_3 + c2_3*v2_3 + c3_3*v3_3;   % P33  

p1_eq = c1_1*v1_1;
p2_eq = c1_2*v1_2;
p3_eq = c1_3*v1_3;

% Additional constraint such that Peq sums to 1.
% eqn10 = c1_1*v1_1 + c1_2*v1_2 + c1_3*v1_3 == 1;
eq_sum = p1_eq + p2_eq + p3_eq == 1;

eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eq_sum];
vars = [c1_1, c1_2, c1_3, c2_1, c2_2, c2_3, c3_1, c3_2, c3_3];

sol = solve(eqns,vars)

p1_eq = (sol.c1_1)*v1_1;
p2_eq = (sol.c1_2)*v1_2;
p3_eq = (sol.c1_3)*v1_3;

c1_1 = sol.c1_1;
c1_2 = sol.c1_2;
c1_3 = sol.c1_3;
c2_1 = sol.c2_1;
c2_2 = sol.c2_2;
c2_3 = sol.c2_3;
c3_1 = sol.c3_1;
c3_2 = sol.c3_2;
c3_3 = sol.c3_3;

save('weightingCoeff_solved_3state123.mat','p1_eq','p2_eq','p3_eq',...
    'c1_1','c1_2','c1_3','c2_1','c2_2','c2_3','c3_1','c3_2','c3_3')

% Define the eigenvector components
syms v0_1 v0_2 v0_3 % compoents of zeroth eigenvector
syms v1_1 v1_2 v1_3 %Components of the first eigenvector
syms v2_1 v2_2 v2_3 %Components of the second eigenvector

% Define expansion coefficients c's for our constants
syms c0_1 c0_2 c0_3  %Components of the zeroth coefficient
syms c1_1 c1_2 c1_3 %Components of the first coefficient
syms c2_1 c2_2 c2_3 %Components of the second coefficient

% Define the eqilibrium probabilities
syms p1_eq p2_eq p3_eq %Equilibrium Populations

% probability from i to j
% pji(t=0) = c0_i*v0_j + c1_i*v1_j + c2_i*v2_j
% P i --> 3                         
eqn1 = 1 == c0_3*v0_3 + c1_3*v1_3 + c2_3*v2_3;   % P33  
eqn2 = 0 == c0_2*v0_3 + c1_2*v1_3 + c2_2*v2_3;   % P32
eqn3 = 0 == c0_1*v0_3 + c1_1*v1_3 + c2_1*v2_3;   % P31

% P i --> 2
eqn4 = 0 == c0_3*v0_2 + c1_3*v1_2 + c2_3*v2_2;   % P23
eqn5 = 1 == c0_2*v0_2 + c1_2*v1_2 + c2_2*v2_2;   % P22
eqn6 = 0 == c0_1*v0_2 + c1_1*v1_2 + c2_1*v2_2;   % P21

% P i --> 1
eqn7 = 0 == c0_3*v0_1 + c1_3*v1_1 + c2_3*v2_1;   % P13
eqn8 = 0 == c0_2*v0_1 + c1_2*v1_1 + c2_2*v2_1;   % P12
eqn9 = 1 == c0_1*v0_1 + c1_1*v1_1 + c2_1*v2_1;   % P11

p1_eq = c0_1*v0_1;
p2_eq = c0_2*v0_2;
p3_eq = c0_3*v0_3;

% Additional constraint such that Peq sums to 1.
% eqn10 = c0_1*v0_1 + c0_2*v0_2 + c0_3*v0_3 == 1;
eq_sum = p1_eq + p2_eq + p3_eq == 1;

eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eq_sum];
vars = [c0_1, c0_2, c0_3, c1_1, c1_2, c1_3, c2_1, c2_2, c2_3];


tic
sol = solve(eqns,vars)
toc

p1_eq = (sol.c0_1)*v0_1;
p2_eq = (sol.c0_2)*v0_2;
p3_eq = (sol.c0_3)*v0_3;

C0_1 = sol.c0_1;
C0_2 = sol.c0_2;
C0_3 = sol.c0_3;
C1_1 = sol.c1_1;
C1_2 = sol.c1_2;
C1_3 = sol.c1_3;
C2_1 = sol.c2_1;
C2_2 = sol.c2_2;
C2_3 = sol.c2_3;

save('weightingCoeff_solved_3state123.mat','p1_eq','p2_eq','p3_eq',...
    'C0_1','C0_2','C0_3','C1_1','C1_2','C1_3','C2_1','C2_2','C2_3')

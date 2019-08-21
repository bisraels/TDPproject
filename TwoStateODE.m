% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: August 2019
%
% PURPOSE: Calculate the two state ODE  for  the master equation
%
% MODIFICATIONS:

clear all
clc

syms k12 k21 P1(t) P2(t)

% Define rate matrix and population vector (as a function of time)
K = [-k12, k21; k12, -k21];
P(t) = [P1(t); P2(t)];

[Vec, lambda] = eig(K);

eval1 = lambda(1,1);
eval2 = lambda(2,2);

evec1 = Vec(:,1);
evec2 = Vec(:,2);

% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);

% Solve the equations and give an output
Psol = dsolve(eqn);

% Rename output as each population
P1(t) = Psol.P1;
P2(t) = Psol.P2;

syms C1 C2

% Now we want to solve using the initial conditions
% Condition 1: P(t) = [P1(t) = 1,P2(t) = 0]

Csols_1 = solve(P1(0) == 1, P2(0)==0,C1,C2);

% Cij, where i = initial condition (in this case 1), j = number of constant
% C11 = Csols_1.C1;
% C12 = Csols_1.C2;
C11 = Csols_1.C2;
C12 = Csols_1.C1;

% Condition 2: P(t) = [0, 1]

Csols_2 = solve(P1(0) == 0, P2(0)==1,C1,C2);

% Cij, where i = initial condition (in this case 2), j = number of constant
% C21 = Csols_2.C1;
% C22 = Csols_2.C2;
C21 = Csols_2.C2;
C22 = Csols_2.C1;


% Write out conditional probabilities
P11(t) = C11 * evec1(1) * exp(eval1 * t) + C12 * evec2(1) * exp(eval2 * t);
P12(t) = C11 * evec1(2) * exp(eval1 * t) + C12 * evec2(2) * exp(eval2 * t);

P21(t) = C21 * evec1(1) * exp(eval1 * t) + C22 * evec2(1) * exp(eval2 * t);
P22(t) = C21 * evec1(2) * exp(eval1 * t) + C22 * evec2(2) * exp(eval2 * t);

save('symCondProb_2state.mat','P11','P12','P21','P22')



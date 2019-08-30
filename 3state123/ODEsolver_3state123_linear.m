% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019
%
% PURPOSE:  Calculate the Linear Three state ODE  for  the master equation
%           1 <-> 2 <-> 3 
%
% OUTPUT: (1) Conditinoal probabilities {Pij(t)}: symCondProb_3state123_linear.mat
%
% MODIFICATIONS: modified from ODEsolver_4state0123.m
%--------------------------------------------------------------------------

clear all
clc
tic

%Declare all the symbolic variables necessary for the ODE Solver
syms k12 k21 k23 k32
syms P1(t) P2(t) P3(t)

% Define the Rate Matrix K
K = [-k12, k21, 0;...
    k12, (-k21 - k23 ), k32;...
    0, k23, -k32;];

% Make a column vector of the equilibrium Populations
P(t) = [P1(t); P2(t);P3(t)];

% Detailed balance condition
%Not existant Here: Only use for loops

% Solve the eigenvalue problem for the eigenvector and eigenvalue
% These will be complex functions of the extant rates in the system
[Vec, lambda] = eig(K);
eval1 = lambda(1,1);        % First Eigenvalue should always be zero.
eval2 = lambda(2,2);
eval3 = lambda(3,3);

evec1 = Vec(:,1);
evec2 = Vec(:,2);
evec3 = Vec(:,3);

% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);
%OUTPUT: eqn =
%  diff(P1(t), t) == k21*P2(t) - k12*P1(t)
%  diff(P2(t), t) == k12*P1(t) - P2(t)*(k21 + k23) + k32*P3(t)
%  diff(P3(t), t) == k23*P2(t) - k32*P3(t)

% Solve the equations and give an output 
%Pij(t) is prob from i--> j, assuming you start in state i: Pi(t=0)=100%=1 
% Condition 1: P(t) = [P1(t) = 1, P2(t) = 0, P3(t) = 0]
Psol_1 = dsolve(eqn,[P1(0) == 1 , P2(0) == 0 , P3(0) == 0]);
%vpa(x) uses variable-precision floating-point arithmetic (VPA)
P11(t) = vpa(Psol_1.P1);
P12(t) = vpa(Psol_1.P2);
P13(t) = vpa(Psol_1.P3);

% Condition 2: P(t) = [P1(t) = 0, P2(t) = 1, P3(t) = 0]
Psol_2 = dsolve(eqn,[P1(0) == 0 , P2(0) == 1 , P3(0) == 0]);
P21(t) = vpa(Psol_2.P1);
P22(t) = vpa(Psol_2.P2);
P23(t) = vpa(Psol_2.P3);

% Condition 3: P(t) = [P1(t) = 0, P2(t) = 0, P3(t) = 1]
Psol_1 = dsolve(eqn,[P1(0) == 0 , P2(0) == 0 , P3(0) == 1]);
P31(t) = vpa(Psol_1.P1);
P32(t) = vpa(Psol_1.P2);
P33(t) = vpa(Psol_1.P3);

save('symCondProb_3state123_linear.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')

elapsedTime = toc;
task_str = 'Calculate and save the eigenvalues and conditional probabilities';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);


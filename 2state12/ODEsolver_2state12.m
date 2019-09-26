% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: August 2019
%
% PURPOSE: Calculate the two state ODE  for  the master equation
%
% INPUT: Nothing. Step 1 of a 2 code process
%
% OUTPUT: (1) Conditional probabilities {Pij(t)}: symCondProb_3state123_linear.mat
%
% MODIFICATIONS:
%
% NOTE: This is updated to use the same method as the 4 state system, but
% the other two data model code does work too. It is easier to understand
% the math from the first two state model code.
%--------------------------------------------------------------------------

clear all;
clc
tic

%--------------------------------------------------------------------------
%  User Prefrences
%--------------------------------------------------------------------------
verbose_mode = 1;
saveMode = 1;

%--------------------------------------------------------------------------
% Solve for the conditional probabilities
%--------------------------------------------------------------------------
%Declare all the symbolic variables necessary for the ODE Solver
syms k12 k21 
syms P1(t) P2(t)

% Define rate matrix and population vector (as a function of time)
K = [-k12 , k21; ...
    k12   , -k21];

if verbose_mode
    disp('Sovling for the eigenvalues of the rate matrix K');
    disp(K);
end

% Make a column vector of the equilibrium Populations
P(t) = [P1(t); P2(t)];

% Solve the eigenvalue problem for the eigenvector and eigenvalue
% These will be complex functions of the extant rates in the system
[Vec, lambda] = eig(K);

eval1 = lambda(1,1);
eval2 = lambda(2,2);

if verbose_mode
    fprintf('Eigenvalue 1 = %s\r',eval1);
    fprintf('Eigenvalue 2 = %s\r',eval2);
end

evec1 = Vec(:,1);
evec2 = Vec(:,2);

if verbose_mode
    fprintf('Eigenvector 1 = (%s,%s)\r',evec1(1),evec1(2));
    fprintf('Eigenvector 2 = (%s,%s)\r',evec2(1),evec2(2));
end

% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);

% Solve the equations and give an output
% Condition 1: when all the population starts in state 1
Psol_1 = dsolve(eqn, P1(0) == 1, P2(0)==0);

% Rename output as conditional probability Pi-->j where i is the initial
% condition,
P11(t) = Psol_1.P1;
P12(t) = Psol_1.P2;


% Condition 2: when all the population starts in state 2
Psol_2 = dsolve(eqn, P1(0) == 0, P2(0)==1);

% Rename output as conditional probability Pi-->j where i is the initial
% condition,
P21(t) = Psol_2.P1;
P22(t) = Psol_2.P2;

%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'calculate the conditional probabilities.';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);

%--------------------------------------------------------------------------
% Save the output
%--------------------------------------------------------------------------
if saveMode
save('symCondProb_2state12.mat','P11','P12','P21','P22')
end


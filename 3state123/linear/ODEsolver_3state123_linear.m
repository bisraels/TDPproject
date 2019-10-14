% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019 (ODEsolver_3state123_linear.m)
%
% PURPOSE:  Calculate the Linear Three state ODE  for  the master equation
%           1 <-> 2 <-> 3 
%
% INPUT: Nothing. Step 1 of a 2 code process
%
% OUTPUT: (1) Conditinoal probabilities {Pij(t)}: symCondProb_3state123_linear.mat
%
% MODIFICATIONS: modified from ODEsolver_4state0123.m
%--------------------------------------------------------------------------

clear all
clc
tic

%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
verbose_mode = 1;


%--------------------------------------------------------------------------
% Solve for the conditional probabilities
%--------------------------------------------------------------------------
%Declare all the symbolic variables necessary for the ODE Solver
syms k12 k21 k23 k32
syms P1(t) P2(t) P3(t)

% Define the Rate Matrix K
K = [-k12, k21, 0;...
    k12, (-k21 - k23 ), k32;...
    0, k23, -k32;];
if verbose_mode
    disp('Sovling for the eigenvalues of the rate matrix K');
    disp(K);
end

% Make a column vector of the equilibrium Populations
P(t) = [P1(t); P2(t); P3(t)];

% Detailed balance condition
%Not existant Here: Only use for loops

% Solve the eigenvalue problem for the eigenvector and eigenvalue
% These will be complex functions of the extant rates in the system
[Vec, lambda] = eig(K);
eval1 = lambda(1,1);        % First Eigenvalue should always be zero.
eval2 = lambda(2,2);
eval3 = lambda(3,3);

if verbose_mode
    fprintf('Eigenvalue 1 = %s\r',eval1);
    fprintf('Eigenvalue 2 = %s\r',eval2);
    fprintf('Eigenvalue 3 = %s\r',eval3);
end

%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'calculate the eigenvalues.';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);


%We do not use these eigenvectors as anything
evec1 = Vec(:,1);
evec2 = Vec(:,2);
evec3 = Vec(:,3);

if verbose_mode
    fprintf('Eigenvector 1 = %s\r',evec1);
    fprintf('Eigenvector 2 = %s\r',evec2);
    fprintf('Eigenvector 3 = %s\r',evec3);
end

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
%The solution to the eigenvector problem will be included by dsove
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

%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'calculate the conditional probabilities.';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);

%--------------------------------------------------------------------------
% Save the output
%--------------------------------------------------------------------------
save('symCondProb_3state123_linear.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')

%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'save the conditional probabilities and eigenvalues.';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);

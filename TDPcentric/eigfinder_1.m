%__________________________________________________________________________
% AUTHOR: Brett Israels
%
% NAME: eigfinder_1.m
%
% FUNCTION: Solves eigenvalues problem of the master equation
%
% PROCEDURE: 
%%%%%%%%%%%%
% [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
%
%%%%%%%
% INPUT:
% (1)
%
% OUTPUT:
% (1) V are the eigenvectors
% (2) D are the eigenvalues
%
% MODIFICATION LOG:
% BI 20190814 Creation
%__________________________________________________________________________


% 1 <--> 2
% dP1/dt = P2eq*k21 - P1eq*k12;
% DP2/dt = P1eqk12 - P2eqk21
%%
syms a b c x
eqn = a*x^2 + b*x + c == 0;
sol = solve(eqn)
sola = solve(eqn, a)
solc = solve(eqn, c)
solx = solve(eqn, x)
% syms k12 k21 P1 P2 
 %%
  syms k12 k21
 M = [-k12, k21;...
     k12, -k21];
[V,D] = eig(M);
disp('The Eigenvalues are: ')
for i = 1:length(D)
    disp(['Eigenvalue ' num2str(i) ' = ']);
disp(D(i,i))
end

disp('The Eigenvectors are: ')
for i = 1:length(V)
    disp(['Eigenvector ' num2str(i) ' = ']);
disp(V(:,i))
end
%%
% 1--2--3--1 cyclical 3 state
syms k12 k13 k21 k23 k32 k31
M = [(-k12-k13) , k21 , k31;...
    k12 , (-k21-k23) , k32;...
    k13 , k23 , (-k31-k32);]

[V,D] = eig(M);
disp('The Eigenvalues are: ')
for i = 1:length(D)
    disp(['Eigenvalue ' num2str(i) ' = ']);
disp(D(i,i))
end

disp('The Eigenvectors are: ')
for i = 1:length(V)
    disp(['Eigenvector ' num2str(i) ' = ']);
disp(V(:,i))
end
%%
% 0-1-2-0 cyclical 3 state
syms k01 k02 k10 k12 k21 k20
M = [(-k01-k02) , k10 , k20;...
    k01 , (-k10-k12) , k21;...
    k02 , k12 , (-k20-k21);]

[V,D] = eig(M);
disp('The Eigenvalues are: ')
for i = 1:length(D)
    disp(['Eigenvalue ' num2str(i) ' = ']);
disp(D(i,i))
end

disp('The Eigenvectors are: ')
for i = 1:length(V)
    disp(['Eigenvector ' num2str(i) ' = ']);
disp(V(:,i))
end
%%
% syms k12 k21 P1 P2 
% eqns = [P1 + P2 == 1, P1*k12 - P2*k21 == 0, P2*k21 - P1*k12 == 0];
% vars = [k12 k21 P1 P2];
% [solk12, solk21, solP1, solP2] = solve(eqns, vars);

%%
syms k12 k02 k10 k20 k01 k12 k21 k23 k32 
M = [-k01-k02,k10,k20,0;k01,-k10-k12,k21,0;k02,k12,-k20-k21-k23,k32;0,0,k23,-k32];
[V,D] = eig(M);





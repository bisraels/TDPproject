% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019
%
% PURPOSE:  Calculate the five state ODE for the master equation (where 2
%           and 2' are different positions of the dimer on the ssDNA)
%           1 <-> 0 <-> 1' <-> 2 <-> 2' <-> 1'
%           Note: model labeling    0   1   1'  2   2'
%                 code labeling     0   1   2   3   4
%           This is adapted from the FourStateODE_0123_fixed.m
%
% MODIFICATIONS:
%--------------------------------------------------------------------------

clear all
clc
tic
% syms k01 k02 k10 k20 k23 k24 k32 k34 k42 k43 
 syms P0(t) P1(t) P2(t) P3(t) P4(t)

t = sym('t');

%--------------------------------------------------------------------------
% Define rates:
%--------------------------------------------------------------------------
k01 = 1;
k02 = 2;
k03 = 3;    % k03 = 0;
k04 = 4;    % k04 = 0;
k10 = 5;
k12 = 6;    % k12 = 0;
k13 = 7;    % k13 = 0;
k14 = 8;    % k14 = 0;
k20 = 9;
k21 = 10;    % k21 = 0;
k23 = 11;
k24 = 12;
k30 = 13;    % k30 = 0;
k31 = 14;    % k31 = 0;
% k32 = 87;
k34 = 15;
k40 = 16;    % k40 = 0;
k41 = 17;    % k41 = 0;
k42 = 18;
k43 = 19;
% NOTE: I put in the rates we know are zero in by hand. For this model we
% have 10 non-zero rates.

% Detailed balance condition
k32 = (k23 * k34 * k42)/(k24 * k43);

% Define rate matrix
K = [-(k01+k02+k03+k04), k10,             k20,              k30,              k40;...
        k01,      -(k10+k12+k13+k14),     k21,              k31,              k41;...
        k02,             k12,       -(k20+k21+k23+k24),     k32,              k42;...
        k03,             k13,             k23,      -(k30+k31+k32+k34),       k43;...
        k04,             k14,             k24,              k34,        -(k40+k41+k42+k43)];
% Vector of populations (or probabilities)
P(t) = [P0(t); P1(t); P2(t); P3(t); P4(t)];

[Vec, lambda] = eig(vpa(K));

[lambdaSort, index] = sort(diag(lambda),'descend');
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);

% eval0 = lambda(1,1);        % eval0 should always be zero.
% eval1 = lambda(2,2);
% eval2 = lambda(3,3);
% eval3 = lambda(4,4);
% eval4 = lambda(5,5);
eval0 = double(lambdaSorted(1,1));        % eval0 should always be zero.
if eval0 == 0
disp('First eigenvalue is zero!')
end
% If matlab says the first eigenvalue is barely not zero, force it to be
% zero.
if abs(eval0)< (1e-20)
    eval0 = 0;
    disp('Setting first eigenvalue to zero, because magnitude is less than 1e-20.')
end
eval1 = double(lambdaSorted(2,2));
eval2 = double(lambdaSorted(3,3));
eval3 = double(lambdaSorted(4,4));
eval4 = double(lambdaSorted(5,5));

% evec0 = Vec(:,1);
% evec1 = Vec(:,2);
% evec2 = Vec(:,3);
% evec3 = Vec(:,4);
% evec4 = Vec(:,5);
evec0 = double(VecSorted(:,1));
evec1 = double(VecSorted(:,2));
evec2 = double(VecSorted(:,3));
evec3 = double(VecSorted(:,4));
evec4 = double(VecSorted(:,5));

toc    % took 87 sec with k's symbolic
       % took 1.6 second with k's doubles, but P's symbolic
%%
tic
% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);

% Solve the equations and give an output
% Condition 4: P(t=0) = [P0(0) = 1, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 0]
Psol_0 = dsolve(eqn,[P0(0)==1 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==0]);
P00(t) = real(vpa(Psol_0.P0));
P01(t) = real(vpa(Psol_0.P1));
P02(t) = real(vpa(Psol_0.P2));
P03(t) = real(vpa(Psol_0.P3));
P04(t) = real(vpa(Psol_0.P4));


% Condition 1: P(t=0) = [P0(0) = 0, P1(0) = 1, P2(0) = 0, P3(0) = 0, P4(0) = 0]
Psol_1 = dsolve(eqn,[P0(0) == 0 P1(0)==1 P2(0)==0 P3(0)==0 P4(0)==0]);
P10(t) = real(vpa(Psol_1.P0));
P11(t) = real(vpa(Psol_1.P1));
P12(t) = real(vpa(Psol_1.P2));
P13(t) = real(vpa(Psol_1.P3));
P14(t) = real(vpa(Psol_1.P4));


% Condition 2: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 1, P3(0) = 0, P4(0) = 0]
Psol_2 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==1 P3(0)==0 P4(0)==0]);
P20(t) = real(vpa(Psol_2.P0));
P21(t) = real(vpa(Psol_2.P1));
P22(t) = real(vpa(Psol_2.P2));
P23(t) = real(vpa(Psol_2.P3));
P24(t) = real(vpa(Psol_2.P4));


% Condition 3: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 1, P4(0) = 0]
Psol_3 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==1 P4(0)==0]);
P30(t) = real(vpa(Psol_3.P0));
P31(t) = real(vpa(Psol_3.P1));
P32(t) = real(vpa(Psol_3.P2));
P33(t) = real(vpa(Psol_3.P3));
P34(t) = real(vpa(Psol_3.P4));

% Condition 4: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 1]
Psol_4 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==1]);
P40(t) = real(vpa(Psol_4.P0));
P41(t) = real(vpa(Psol_4.P1));
P42(t) = real(vpa(Psol_4.P2));
P43(t) = real(vpa(Psol_4.P3));
P44(t) = real(vpa(Psol_4.P4));

% Do the equilibrium probabilities sum  to 1?
if P00(inf)+P11(inf)+P22(inf)+P33(inf)+P44(inf) == 1
    disp('Equlibrium probabilities sum to 1!')
else
    disp('Equilibrium probabilities do not sum to 1.')
    X = ['Instead they sum to %d', P00(inf)+P11(inf)+P22(inf)+P33(inf)+P44(inf)];
    disp(X)
end


save('symCondProb_5state01234.mat','P00','P01','P02','P03','P04','P10','P11','P12','P13','P14','P20','P21','P22','P23','P24','P30','P31','P32','P33','P34','P40','P41','P42','P43','P44','eval0','eval1','eval2','eval3','eval4')

toc



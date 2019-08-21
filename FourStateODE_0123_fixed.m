% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019
%
% PURPOSE:  Calculate the four state ODE  for  the master equation
%           0 <-> 1 <-> 1' <-> 0 or 0 <-> 1 <-> 2
%           Note: model labeling    0   1   1'  2
%                 code labeling     0   1   2   3
%           This is the updated version that eliminates the steps that were
%           previously giving us complex probabilities
%
% MODIFICATIONS:
%--------------------------------------------------------------------------

clear all
clc

syms k01 k02 k03 k10 k12 k13 k20 k21 k23 k30 k31 k32 positive 
syms P0(t) P1(t) P2(t) P3(t)
% assume(P0(t),'real')
% assume(P1(t),'real')
% assume(P2(t),'real')
% assume(P3(t),'real')

% Define rate matrix and population vector (as a function of time)
% K = [-(k01 + k02 + k03), k10, k20, k30;...
%     k01, -(k10 + k12 + k13), k21, k31;...
%     k02, k12, -(k20 + k21 + k23), k13;...
%     k03, k13, k23, -(k30 + k31 + k32);];
K = [-(k01 + k02 + 0), k10, k20, 0;...
    k01, -(k10 + k12 + 0), k21, 0;...
    k02, k12, -(k20 + k21 + k23), k32;...
    0, 0, k23, -(0 + 0 + k32);];
P(t) = [P0(t);P1(t); P2(t);P3(t)];

% Detailed balance condition
%k20 = (k10 * k02 * k21)/(k01 * k12);
k02 = (k01 * k12 * k20)/(k10 * k21);

[Vec, lambda] = eig(K);
eval0 = real(lambda(1,1));        % eval0 should always be zero.
eval1 = real(lambda(2,2));
eval2 = real(lambda(3,3));
eval3 = real(lambda(4,4));

evec0 = Vec(:,1);
evec1 = Vec(:,2);
evec2 = Vec(:,3);
evec3 = Vec(:,4);


% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);

% Solve the equations and give an output
% Condition 0: P(t) = [P0(t) = 1, P1(t) = 0, P2(t) = 0, P3(t) = 0]
Psol_0 = dsolve(eqn,[P0(0) == 1 P1(0)==0 P2(0)==0 P3(0)==0]);
P00(t) = vpa(Psol_0.P0);
P01(t) = vpa(Psol_0.P1);
P02(t) = vpa(Psol_0.P2);
P03(t) = vpa(Psol_0.P3);

% Condition 1: P(t) = [P0(t) = 0, P1(t) = 1, P2(t) = 0, P3(t) = 0]
Psol_1 = dsolve(eqn,[P0(0) == 0 P1(0)==1 P2(0)==0 P3(0)==0]);
P10(t) = vpa(Psol_1.P0);
P11(t) = vpa(Psol_1.P1);
P12(t) = vpa(Psol_1.P2);
P13(t) = vpa(Psol_1.P3);

% Condition 2: P(t) = [P0(t) = 0, P1(t) = 0, P2(t) = 1, P3(t) = 0]
Psol_2 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==1 P3(0)==0]);
P20(t) = vpa(Psol_2.P0);
P21(t) = vpa(Psol_2.P1);
P22(t) = vpa(Psol_2.P2);
P23(t) = vpa(Psol_2.P3);

% Condition 3: P(t) = [P0(t) = 0, P1(t) = 0, P2(t) = 0, P3(t) = 3]
Psol_3 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==1]);
P30(t) = vpa(Psol_3.P0);
P31(t) = vpa(Psol_3.P1);
P32(t) = vpa(Psol_3.P2);
P33(t) = vpa(Psol_3.P3);

save('symCondProb_4state0123_fixed.mat','P00','P01','P02','P03','P10','P11','P12','P13','P20','P21','P22','P23','P30','P31','P32','P33')



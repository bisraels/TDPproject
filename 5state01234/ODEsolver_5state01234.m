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
syms k01 k02 k10 k20 k23 k24 k32 k34 k42 k43 
syms P0(t) P1(t) P2(t) P3(t) P4(t)

K = [-(k01 + k02),      k10,            k20,        0,              0;
        k01,          -(k10),            0,         0,              0;
        k02,            0,       -(k20+k23+k24),   k32,            k42;
        0,              0,              k23,    -(k32+k42),         k43;
        0,              0               k24,        k34,        -(k42+k43)];
P(t) = [P0(t);P1(t); P2(t);P3(t); P4(t)];

% Detailed balance condition
k32 = (k23 * k34 * k42)/(k24 * k43);

[Vec, lambda] = eig(K);
eval0 = lambda(1,1);        % eval0 should always be zero.
eval1 = lambda(2,2);
eval2 = lambda(3,3);
eval3 = lambda(4,4);
eval4 = lambda(5,5);

evec0 = Vec(:,1);
evec1 = Vec(:,2);
evec2 = Vec(:,3);
evec3 = Vec(:,4);
evec4 = Vec(:,5);
toc    % took 87 sec
%%
tic
% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);

% Solve the equations and give an output
% Condition 4: P(t=0) = [P0(0) = 1, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 0]
tic
Psol_0 = dsolve(eqn,[P0(0) == 1 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==0]);
P00(t) = vpa(Psol_0.P0);
P01(t) = vpa(Psol_0.P1);
P02(t) = vpa(Psol_0.P2);
P03(t) = vpa(Psol_0.P3);
P04(t) = vpa(Psol_0.P4);
toc

% Condition 1: P(t=0) = [P0(0) = 0, P1(0) = 1, P2(0) = 0, P3(0) = 0, P4(0) = 0]
Psol_1 = dsolve(eqn,[P0(0) == 0 P1(0)==1 P2(0)==0 P3(0)==0 P4(0)==0]);
P10(t) = vpa(Psol_1.P0);
P11(t) = vpa(Psol_1.P1);
P12(t) = vpa(Psol_1.P2);
P13(t) = vpa(Psol_1.P3);
P14(t) = vpa(Psol_1.P4);


% Condition 2: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 1, P3(0) = 0, P4(0) = 0]
Psol_2 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==1 P3(0)==0 P4(0)==0]);
P20(t) = vpa(Psol_2.P0);
P21(t) = vpa(Psol_2.P1);
P22(t) = vpa(Psol_2.P2);
P23(t) = vpa(Psol_2.P3);
P24(t) = vpa(Psol_2.P4);


% Condition 3: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 1, P4(0) = 0]
Psol_3 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==1 P4(0)==0]);
P30(t) = vpa(Psol_3.P0);
P31(t) = vpa(Psol_3.P1);
P32(t) = vpa(Psol_3.P2);
P33(t) = vpa(Psol_3.P3);
P34(t) = vpa(Psol_3.P4);

% Condition 4: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 1]
Psol_3 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==1]);
P40(t) = vpa(Psol_4.P0);
P41(t) = vpa(Psol_4.P1);
P42(t) = vpa(Psol_4.P2);
P43(t) = vpa(Psol_4.P3);
P44(t) = vpa(Psol_4.P4);


save('symCondProb_5state01234_fixed.mat','P00','P01','P02','P03','P04','P10','P11','P12','P13','P14''P20','P21','P22','P23','P24','P30','P31','P32','P33','P34','P40','P41','P42','P43','P44','eval0','eval1','eval2','eval3','eval4')

toc



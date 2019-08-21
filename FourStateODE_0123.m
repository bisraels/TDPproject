% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: August 2019
%
% PURPOSE: Calculate the four state ODE  for  the master equation
%           0 <-> 1 <-> 1' <-> 0 or 0 <-> 1 <-> 2
%           Note: model labeling    0   1   1'  2
%                 code labeling     0   1   2   3
%
% MODIFICATIONS:
% NOTE: THIS CODE HAS PROBLEMS WHERE IT GIVES IMAGINARY PROBABILITIES. USE
% THIS ONE INSTEAD: FourStateODE_0123_fixed

clear all
clc

syms k01 k02 k03 k10 k12 k13 k20 k21 k23 k30 k31 k32 positive 
syms P0(t) P1(t) P2(t) P3(t)
assume(P0(t),'real')
assume(P1(t),'real')
assume(P2(t),'real')
assume(P3(t),'real')

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

%%
% Try looking at left and right evec's
[VecR, lambda, VecL] = eig(K);
%[V, D, W] = eig(K);
eval0 = lambda(1,1);     % eval0 should always be zero.
eval1 = lambda(2,2);
eval2 = lambda(3,3);
eval3 = lambda(4,4);

evecR0 = VecR(:,1);
evecR1 = VecR(:,2);
evecR2 = VecR(:,3);
evecR3 = VecR(:,4);

evecL0 = VecL(:,1);
evecL1 = VecL(:,2);
evecL2 = VecL(:,3);      % Use isequaln(a,b) to compare symbolic expressions
evecL3 = VecL(:,4);      % This tells us that the left and right evec's are the same.

%%
% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);

% Solve the equations and give an output
Psol = dsolve(eqn);

% Rename output as each population
P0(t) = vpa(Psol.P0);   % vpa() forces the computer to find the solution instead of leaving it in terms of roots()
P1(t) = vpa(Psol.P1);
P2(t) = vpa(Psol.P2);
P3(t) = vpa(Psol.P3);


syms C1 C2 C3 C4

% Now we want to solve using the initial conditions
% Condition 0: P(t) = [P0(t) = 1, P1(t) = 0, P2(t) = 0, P3(t) = 0]
Csols_0 = solve(P0(0) == 1, P1(0)==0, P2(0)==0, P3(0)==0, C1, C2, C3, C4);

% Cij, where i = initial condition (in this case state zero starts with 1, 
% so lets call it condition 0), j = number of constant
C00 = vpa(Csols_0.C1);
C01 = vpa(Csols_0.C2);
C02 = vpa(Csols_0.C3);
C03 = vpa(Csols_0.C4);
% C00 = real(vpa(Csols_0.C1));
% C01 = real(vpa(Csols_0.C2));
% C02 = real(vpa(Csols_0.C3));
% C03 = real(vpa(Csols_0.C4));

% Condition 1: P(t) = [P0(t=0) = 0, P1(t=0) = 1, P2(t=0) = 0, P3(t=0) = 0]
Csols_1 = solve(P0(0) == 0, P1(0)== 1, P2(0)== 0, P3(0)== 0, C1, C2, C3, C4);
C10 = vpa(Csols_1.C1);
C11 = vpa(Csols_1.C2);
C12 = vpa(Csols_1.C3);
C13 = vpa(Csols_1.C4);
% C10 = real(vpa(Csols_1.C1));
% C11 = real(vpa(Csols_1.C2));
% C12 = real(vpa(Csols_1.C3));
% C13 = real(vpa(Csols_1.C4));

% Condition 2: P(t) = [P0(t=0) = 0, P1(t=0) = 0, P2(t=0) = 1, P3(t=0) = 0]
Csols_2 = solve(P0(0) == 0, P1(0)== 0, P2(0)== 1, P3(0)== 0, C1, C2, C3, C4);
C20 = vpa(Csols_2.C1);
C21 = vpa(Csols_2.C2);
C22 = vpa(Csols_2.C3);
C23 = vpa(Csols_2.C4);
% C20 = real(vpa(Csols_2.C1));
% C21 = real(vpa(Csols_2.C2));
% C22 = real(vpa(Csols_2.C3));
% C23 = real(vpa(Csols_2.C4));

% Condition 3: P(t) = [P0(t=0) = 0, P1(t=0) = 0, P2(t=0) = 0, P3(t=0) = 1]
Csols_3 = solve(P0(0) == 0, P1(0)== 0, P2(0)== 0, P3(0)== 1, C1, C2, C3, C4);
C30 = vpa(Csols_3.C1);
C31 = vpa(Csols_3.C2);
C32 = vpa(Csols_3.C3);
C33 = vpa(Csols_3.C4);
% C30 = real(vpa(Csols_3.C1));
% C31 = real(vpa(Csols_3.C2));
% C32 = real(vpa(Csols_3.C3));
% C33 = real(vpa(Csols_3.C4));

% Write out conditional probabilities ( labeled as: Pi-->j and i labels the initial condition)
P00(t) = C00 * evec0(1)*exp(eval0 * t) + C01 * evec1(1) * exp(eval1 * t) + C02 * evec2(1) * exp(eval2 * t) + C03 * evec3(1) * exp(eval3 * t);
P01(t) = C00 * evec0(2)*exp(eval0 * t) + C01 * evec1(2) * exp(eval1 * t) + C02 * evec2(2) * exp(eval2 * t) + C03 * evec3(2) * exp(eval3 * t);
P02(t) = C00 * evec0(3)*exp(eval0 * t) + C01 * evec1(3) * exp(eval1 * t) + C02 * evec2(3) * exp(eval2 * t) + C03 * evec3(3) * exp(eval3 * t);
P03(t) = C00 * evec0(4)*exp(eval0 * t) + C01 * evec1(4) * exp(eval1 * t) + C02 * evec2(4) * exp(eval2 * t) + C03 * evec3(4) * exp(eval3 * t);

P10(t) = C10 * evec0(1)*exp(eval0 * t) + C11 * evec1(1) * exp(eval1 * t) + C12 * evec2(1) * exp(eval2 * t) + C13 * evec3(1) * exp(eval3 * t);
P11(t) = C10 * evec0(2)*exp(eval0 * t) + C11 * evec1(2) * exp(eval1 * t) + C12 * evec2(2) * exp(eval2 * t) + C13 * evec3(2) * exp(eval3 * t);
P12(t) = C10 * evec0(3)*exp(eval0 * t) + C11 * evec1(3) * exp(eval1 * t) + C12 * evec2(3) * exp(eval2 * t) + C13 * evec3(3) * exp(eval3 * t);
P13(t) = C10 * evec0(4)*exp(eval0 * t) + C11 * evec1(4) * exp(eval1 * t) + C12 * evec2(4) * exp(eval2 * t) + C13 * evec3(4) * exp(eval3 * t);

P20(t) = C20 * evec0(1)*exp(eval0 * t) + C21 * evec1(1) * exp(eval1 * t) + C22 * evec2(1) * exp(eval2 * t) + C23 * evec3(1) * exp(eval3 * t);
P21(t) = C20 * evec0(2)*exp(eval0 * t) + C21 * evec1(2) * exp(eval1 * t) + C22 * evec2(2) * exp(eval2 * t) + C23 * evec3(2) * exp(eval3 * t);
P22(t) = C20 * evec0(3)*exp(eval0 * t) + C21 * evec1(3) * exp(eval1 * t) + C22 * evec2(3) * exp(eval2 * t) + C23 * evec3(3) * exp(eval3 * t);
P23(t) = C20 * evec0(4)*exp(eval0 * t) + C21 * evec1(4) * exp(eval1 * t) + C22 * evec2(4) * exp(eval2 * t) + C23 * evec3(4) * exp(eval3 * t);

P30(t) = C30 * evec0(1)*exp(eval0 * t) + C31 * evec1(1) * exp(eval1 * t) + C32 * evec2(1) * exp(eval2 * t) + C33 * evec3(1) * exp(eval3 * t);
P31(t) = C30 * evec0(2)*exp(eval0 * t) + C31 * evec1(2) * exp(eval1 * t) + C32 * evec2(2) * exp(eval2 * t) + C33 * evec3(2) * exp(eval3 * t);
P32(t) = C30 * evec0(3)*exp(eval0 * t) + C31 * evec1(3) * exp(eval1 * t) + C32 * evec2(3) * exp(eval2 * t) + C33 * evec3(3) * exp(eval3 * t);
P33(t) = C30 * evec0(4)*exp(eval0 * t) + C31 * evec1(4) * exp(eval1 * t) + C32 * evec2(4) * exp(eval2 * t) + C33 * evec3(4) * exp(eval3 * t);


save('symCondProb_4state0123.mat','P00','P01','P02','P03','P10','P11','P12','P13','P20','P21','P23','P30','P31','P32','P33')



% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019
%
% PURPOSE:  Evaluate the five state (01234) conditional probabilties with a set of 
%           rates then the  2 point TCF and 4 point TCF for that system
%
% METHOD:   Use ode45 to solve for the conditional probabilities
%               - this is different than the other program in that ode45
%               calculates the probabilities over a specific range of
%               times, so we wouldn't need to then define a matrix of
%               vectors of probability over time, this would create each
%               probability as a vector of times.
%               - The down side, is you will need a different way to
%               calculate the equilibrium probabilities, because the
%               program can't calculate up to infinite time. 
%
% INPUT:    Use function: ODE_01234.m
%           
%
% MODIFICATIONS:
%
%__________________________________________________________________________

clear all

t = sym('t');

k01 = 0.001;
k02 = 2;
k03 = 0;
k04 = 0;
k10 = 0.000034;
k12 = 0;
k13 = 0;
k14 = 0;
k20 = 0.153;
k21 = 0;
k23 = 700;
k24 = 34;
k30 = 0;
k31 = 0;
% k32 = 87;
k34 = 45;
k40 = 0;
k41 = 0;
k42 = 4;
k43 = 7;
% Detailed balance condition
k32 = (k23 * k34 * k42)/(k24 * k43);


K = [-(k01+k02+k03+k04), k10,             k20,              k30,              k40;...
        k01,      -(k10+k12+k13+k14),     k21,              k31,              k41;...
        k02,             k12,       -(k20+k21+k23+k24),     k32,              k42;...
        k03,             k13,             k23,      -(k30+k31+k32+k34),       k43;...
        k04,             k14,             k24,              k34,        -(k40+k41+k42+k43)];
    
    syms P0(t) P1(t) P2(t) P3(t) P4(t)

% Vector of populations (or probabilities)
P(t) = [P0(t); P1(t); P2(t); P3(t); P4(t)];


% ODE to solve
% eqn = diff(P(t),t)== K * P(t);

t0 = 599;
tf = 600;
tspan = [t0 tf];
tic

cond0 = [1 0 0 0 0];
cond1 = [0 1 0 0 0];
cond2 = [0 0 1 0 0];
cond3 = [0 0 0 1 0];
cond4 = [0 0 0 0 1];
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

[~,Psol_0] = ode45('ODE_01234',tspan,cond0,opts,k01,k02,k03,k04,k10,k12,k13,k14,k20,k21,k23,k24,k30,k31,k32,k34,k40,k41,k42,k43);
P00 = Psol_0(:,1);
P01 = Psol_0(:,2);
P02 = Psol_0(:,3);
P03 = Psol_0(:,4);
P04 = Psol_0(:,5);

[~,Psol_1] = ode45('ODE_01234',tspan,cond1,opts,k01,k02,k03,k04,k10,k12,k13,k14,k20,k21,k23,k24,k30,k31,k32,k34,k40,k41,k42,k43);
P10 = Psol_1(:,1);
P11 = Psol_1(:,2);
P12 = Psol_1(:,3);
P13 = Psol_1(:,4);
P14 = Psol_1(:,5);

[~,Psol_2] = ode45('ODE_01234',tspan,cond2,opts,k01,k02,k03,k04,k10,k12,k13,k14,k20,k21,k23,k24,k30,k31,k32,k34,k40,k41,k42,k43);
P20 = Psol_2(:,1);
P21 = Psol_2(:,2);
P22 = Psol_2(:,3);
P23 = Psol_2(:,4);
P24 = Psol_2(:,5);

[~,Psol_3] = ode45('ODE_01234',tspan,cond3,opts,k01,k02,k03,k04,k10,k12,k13,k14,k20,k21,k23,k24,k30,k31,k32,k34,k40,k41,k42,k43);
P30 = Psol_3(:,1);
P31 = Psol_3(:,2);
P32 = Psol_3(:,3);
P33 = Psol_3(:,4);
P34 = Psol_3(:,5);

[~,Psol_4] = ode45('ODE_01234',tspan,cond4,opts,k01,k02,k03,k04,k10,k12,k13,k14,k20,k21,k23,k24,k30,k31,k32,k34,k40,k41,k42,k43);
P40 = Psol_4(:,1);
P41 = Psol_4(:,2);
P42 = Psol_4(:,3);
P43 = Psol_4(:,4);
P44 = Psol_4(:,5);


% Maybe equilibrium populations:
    disp('End values sum to')
    disp(P00(length(P00)) + P11(length(P11)) + P22(length(P22)) + P33(length(P33)) + P44(length(P44)) )

    

toc




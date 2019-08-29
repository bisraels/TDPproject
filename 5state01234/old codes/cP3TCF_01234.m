% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019
%
% PURPOSE:  Evaluate the five state (01234) conditional probabilties with a set of 
%           rates then the  2 point TCF and 4 point TCF for that system
%
% METHOD:    Define rates, detailed balance, & evaluate eigenvalues
%               (1) Conditional probabilites, Pji 
%               (2) Equilibrium probabilites, Pi_eq 
%               (3) FRET states, Ai 
%            Calculate:
%               (1) Two point TCF, C2
%               (2) Four point TCF, C4
%
% INPUT:    conditional probabilities from ODEsolver_5state01234.m 
%           loaded from: 'symCondProb_5state01234.mat'
%           
%
% MODIFICATIONS:
%
%__________________________________________________________________________

%--------------------------------------------------------------------------
% LOAD VARIABLES
%--------------------------------------------------------------------------
% Load output created by ODEsolver_5state01234.m
disp('Loading the conditional Probabilities as a function of rates');
tic
%load('symCondProb_4state0123.mat')        % This is from the old code
load('symCondProb_5state01234.mat')   % This is the output of the fixed code.
toc

tic
disp('Calculating conditional probabilities using the rates defined')
t = sym('t');

%--------------------------------------------------------------------------
% Define rates:
%--------------------------------------------------------------------------
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
% NOTE: I put in the rates we know are zero in by hand. For this model we
% have 10 non-zero rates.

% Detailed balance condition
k32 = (k23 * k34 * k42)/(k24 * k43);

%--------------------------------------------------------------------------
% EIGENVALUES
%--------------------------------------------------------------------------
% Evaluate the eigenvalues in terms of the rates defined above - produce as
% doubles
eval0 = double(vpa(subs(eval0)));
eval1 = double(vpa(subs(eval1)));
eval2 = double(vpa(subs(eval2)));
eval3 = double(vpa(subs(eval3)));
eval4 = double(vpa(subs(eval4)));

%--------------------------------------------------------------------------
% CONDITIONAL PROBABILITIES
%--------------------------------------------------------------------------
% Evaluate conditional probabilties by substituting in values from above
% and using vpa() to force the simplest form of the output.

t = sym('t');   % For now, the conditional probabilities are a function of t.

P00(t) = vpa(subs(P00));    
P01(t) = vpa(subs(P01));
P02(t) = vpa(subs(P02));
P03(t) = vpa(subs(P03));

P10(t) = vpa(subs(P10));
P11(t) = vpa(subs(P11));
P12(t) = vpa(subs(P12));
P13(t) = vpa(subs(P13));

P20(t) = vpa(subs(P20));
P21(t) = vpa(subs(P21));
P22(t) = vpa(subs(P22));
P23(t) = vpa(subs(P23));

P30(t) = vpa(subs(P30));
P31(t) = vpa(subs(P31));
P32(t) = vpa(subs(P32));
P33(t) = vpa(subs(P33));
% Note: these output as symbolic functions ('symfun') need to evaluate at a
% time and call it with double() to turn it into a number.


% Define a matrix of the conditional probabilities: (Pi-->j with i is initial condition)
 cP = [P00(t), P10(t), P20(t), P30(t), P40(t);
       P01(t), P11(t), P21(t), P31(t), P41(t);
       P02(t), P12(t), P22(t), P32(t), P42(t);
       P03(t), P13(t), P23(t), P33(t), P43(t);
       P04(t), P14(t), P24(t), P34(t), P44(t)];
% Row = final state & Column = initial condition

%--------------------------------------------------------------------------
% EQUILIBRIUM PROBABILITIES
%--------------------------------------------------------------------------
% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P0EQ = P00(inf);
P1EQ = P11(inf);
P2EQ = P22(inf);
P3EQ = P33(inf);
P4EQ = P44(inf);
 
Peq = [P0EQ; P1EQ; P2EQ; P3EQ; P4EQ];

if sum(Peq) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('Problem: Equilibrium probabilities DO NOT sum to 1.')
end

%--------------------------------------------------------------------------
% FRET STATES
%--------------------------------------------------------------------------
% Define FRET states in a vector A:
A0 = 0.34;  % These are just guesses for now
A1 = 0.42;
A2 = 0.51;
A3 = 0.67;
A5 = 0.49;

A = [A0; A1; A2; A3; A4];

% Adapt A so that it will go to the right limits for C2: (difference
% correlation function)
Amean = sum(A.*Peq);
A = A - Amean;

msq = sum((A.^2).*Peq);     % square of mean <A^2>
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2


%% CALCULATE TCF's 

%--------------------------------------------------------------------------
% Two point TCF:
%--------------------------------------------------------------------------

C2(t) = 0*t;
for i = 1:numel(A)
    for j = 1:numel(A)
        
        C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
        C2(t) = C2(t) + C2temp(t);
    end
end


if double(C2(0)) == double(msq)
    disp('Mean of the square <A^2> matches C2(t=0)!')
else 
    disp('Problem: mean of the square <A^2> DOES NOT match C2(t=0)')
end

if double(C2(10^20)) == double(sqm)
    disp('Square of the mean <A>^2 matches C2(t=inf)!')
else 
    disp('Problem: square of the mean <A>^2 DOES NOT match C2(t=inf)')
end

%--------------------------------------------------------------------------
% Plot two point TCF
%--------------------------------------------------------------------------

close all
figure
TCF2pt = fplot(C2(t),[0,1]);
title('Analytical Two point TCF','FontSize',18)
xlabel('Time (\tau_1)','FontSize',14);
ylabel('C^{(2)}(\tau) (\tau)','FontSize',14);
ax = gca;
ax.XScale = 'log';



%%
%--------------------------------------------------------------------------
% Four point TCF:
%--------------------------------------------------------------------------

% Define times:
Npts = 50;
tau1vec = logspace(0,3,Npts);   % Vector of times
tau2 = 1;                       % Tau2 is a constant
tau3vec = tau1vec;              % Tau3 is the same vector of time as tau1

t1 = tau1vec;
t2 = tau2;
t3 = tau3;
% These are what get inserted into the conditional probability functions.


tic
disp('... Calculating the Conditional Probabilities');
% Initialize size of 3D array for conditional probability for t1:
cP_t1 = zeros(5,5,length(t1));

% Define the vector in each position
cP_t1(1,1,:) = double(P00(t1));
cP_t1(2,1,:) = double(P01(t1));
cP_t1(3,1,:) = double(P02(t1));
cP_t1(4,1,:) = double(P03(t1));
cP_t1(5,1,:) = double(P04(t1));

cP_t1(1,2,:) = double(P10(t1));
cP_t1(2,2,:) = double(P11(t1));
cP_t1(3,2,:) = double(P12(t1));
cP_t1(4,2,:) = double(P13(t1));
cP_t1(5,2,:) = double(P14(t1));

cP_t1(1,3,:) = double(P20(t1));
cP_t1(2,3,:) = double(P21(t1));
cP_t1(3,3,:) = double(P22(t1));
cP_t1(4,3,:) = double(P23(t1));
cP_t1(5,3,:) = double(P24(t1));

cP_t1(1,4,:) = double(P30(t1));
cP_t1(2,4,:) = double(P31(t1));
cP_t1(3,4,:) = double(P32(t1));
cP_t1(4,4,:) = double(P33(t1));
cP_t1(5,4,:) = double(P34(t1));

cP_t1(1,5,:) = double(P40(t1));
cP_t1(2,5,:) = double(P41(t1));
cP_t1(3,5,:) = double(P42(t1));
cP_t1(4,5,:) = double(P43(t1));
cP_t1(5,5,:) = double(P44(t1));

% Initialize size of 3D array for conditional probability for t2:
cP_t2 = zeros(4,4,length(t2));
% Define the vector in each position
cP_t2(1,1,:) = double(P00(t2));
cP_t2(2,1,:) = double(P01(t2));
cP_t2(3,1,:) = double(P02(t2));
cP_t2(4,1,:) = double(P03(t2));
cP_t2(5,1,:) = double(P04(t2));

cP_t2(1,2,:) = double(P10(t2));
cP_t2(2,2,:) = double(P11(t2));
cP_t2(3,2,:) = double(P12(t2));
cP_t2(4,2,:) = double(P13(t2));
cP_t2(5,2,:) = double(P14(t2));

cP_t2(1,3,:) = double(P20(t2));
cP_t2(2,3,:) = double(P21(t2));
cP_t2(3,3,:) = double(P22(t2));
cP_t2(4,3,:) = double(P23(t2));
cP_t2(5,3,:) = double(P24(t2));

cP_t2(1,4,:) = double(P30(t2));
cP_t2(2,4,:) = double(P31(t2));
cP_t2(3,4,:) = double(P32(t2));
cP_t2(4,4,:) = double(P33(t2));
cP_t2(5,4,:) = double(P34(t2));

cP_t2(1,5,:) = double(P40(t2));
cP_t2(2,5,:) = double(P41(t2));
cP_t2(3,5,:) = double(P42(t2));
cP_t2(4,5,:) = double(P43(t2));
cP_t2(5,5,:) = double(P44(t2));

% Initialize size of 3D array for conditional probability for t3:
cP_t3 = zeros(5,5,length(t3));
% Define the vector in each position
cP_t3(1,1,:) = double(P00(t3));
cP_t3(2,1,:) = double(P01(t3));
cP_t3(3,1,:) = double(P02(t3));
cP_t3(4,1,:) = double(P03(t3));
cP_t3(5,1,:) = double(P04(t3));

cP_t3(1,2,:) = double(P10(t3));
cP_t3(2,2,:) = double(P11(t3));
cP_t3(3,2,:) = double(P12(t3));
cP_t3(4,2,:) = double(P13(t3));
cP_t3(5,2,:) = double(P14(t3));

cP_t3(1,3,:) = double(P20(t3));
cP_t3(2,3,:) = double(P21(t3));
cP_t3(3,3,:) = double(P22(t3));
cP_t3(4,3,:) = double(P23(t3));
cP_t3(5,3,:) = double(P24(t3));

cP_t3(1,4,:) = double(P30(t3));
cP_t3(2,4,:) = double(P31(t3));
cP_t3(3,4,:) = double(P32(t3));
cP_t3(4,4,:) = double(P33(t3));
cP_t3(5,4,:) = double(P34(t3));

cP_t3(1,5,:) = double(P40(t3));
cP_t3(2,5,:) = double(P41(t3));
cP_t3(3,5,:) = double(P42(t3));
cP_t3(4,5,:) = double(P43(t3));
cP_t3(5,5,:) = double(P44(t3));

disp('Time to calculate Conditional Probabilities:');
toc


disp('... Calculating the 4 point TCF');
tic
C4vec = [];

 %-------------------------------------------------------------------------
 % Iterate over all the Combinations of FRET States
 %-------------------------------------------------------------------------
     for i = 1:numel(A)
         for j = 1:numel(A)
             for k = 1:numel(A)
                 for l = 1:numel(A)
                     
                     C4term_val =  A(l) *squeeze(cP_t3(l,k,:)) * A(k) * cP_t2(k,j) * A(j) * squeeze(cP_t1(j,i,:))'* A(i) * Peq(i);
                  C4mat = C4mat + C4term_val;

                 end
             end
         end
     end
 C4 = C4mat;
disp('Time to calculate the four point TCF (C4):');
toc

%--------------------------------------------------------------------------
% Plot surface of C4
%--------------------------------------------------------------------------
figure
TCF4point = surf(tau1vec, tau3vec, C4);
title('Analytical Four-point TCF: C^{(4)}','FontSize',18)
xlabel('Time (\tau_1)','FontSize',14);
ylabel('Time (\tau_3)','FontSize',14);
zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';


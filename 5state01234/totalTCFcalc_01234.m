% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019
%
% PURPOSE:  Evaluate the five state (01234) conditional probabilties with a set of 
%           rates then the  2 point TCF and 4 point TCF for that system
%
% METHOD:       (1) Define rates and detailed balance
%               (2) Solve eigenvalue problem
%               (3) Solve rate equation ODE
%               (4) Calculate conditional probabilites, Pji 
%               (5) Equilibrium probabilites, Pi_eq 
%               (6) FRET states, Ai 
%            Calculate:
%               (1) Two point TCF, C2
%               (2) Four point TCF, C4
%
% INPUT:    nothing for now, it is completely self contained.
%           
%
% MODIFICATIONS:
%
%__________________________________________________________________________

%--------------------------------------------------------------------------
% LOAD VARIABLES
%--------------------------------------------------------------------------
% Load output created by ODEsolver_5state01234.m
% disp('Loading the conditional Probabilities as a function of rates');
% tic
% %load('symCondProb_4state0123.mat')        % This is from the old code
% load('symCondProb_5state01234.mat')   % This is the output of the fixed code.
% toc

clear all       % CAREFUL: sometimes if the variables are not cleared the eigenvalue 
clc             % calculated from eig() and while solving the ODE will not match.
                % Then all of the subsitutions will break, this is probably
                % what happened if sum(Peq) = 0.

tic
t = sym('t');
syms P0(t) P1(t) P2(t) P3(t) P4(t)

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
% EIGENVALUE PROBLEM
%--------------------------------------------------------------------------
tic
disp('...Solving eigenvalue problem')
% Define rate matrix
K = [-(k01+k02+k03+k04), k10,             k20,              k30,              k40;...
        k01,      -(k10+k12+k13+k14),     k21,              k31,              k41;...
        k02,             k12,       -(k20+k21+k23+k24),     k32,              k42;...
        k03,             k13,             k23,      -(k30+k31+k32+k34),       k43;...
        k04,             k14,             k24,              k34,        -(k40+k41+k42+k43)];
% Vector of populations (or probabilities)
P(t) = [P0(t); P1(t); P2(t); P3(t); P4(t)];


[Vec, lambda] = eig(vpa(K,10));     % the vpa() should make the vec and lambda output as doubles.

[lambdaSort, index] = sort(diag(lambda),'descend');   % sort just in case, it the program seems to output, closest to 0 -> most negative  
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);

%--------------------------------------------------------------------------
% EIGENVALUES
%--------------------------------------------------------------------------
eval0 = double(vpa(lambdaSorted(1,1),9));        % eval0 should always be zero.
if eval0 == 0
disp('First eigenvalue is zero!')
else 
    % If first eigenvalue is not zero, it will screw up the rest of the
    % probability expressions.
    disp('PROBLEM: need to fix first eigenvalue AND probability expressions')
    disp('eval0 =')
    disp(eval0)
end
eval1 = double(lambdaSorted(2,2));
eval2 = double(lambdaSorted(3,3));
eval3 = double(lambdaSorted(4,4));
eval4 = double(lambdaSorted(5,5));

if eval0<=0 && eval1<=0 && eval2 <=0 && eval3 <=0 && eval4 <=0      % For the exponentials to decay, need all eigenvalues to be < 0.
    disp('All eigenvalues are less than zero!')                     % If this is not true, there is a problem.
else
    disp('PROBLEM: not all eigenvalues are less than zero.')
    disp('Instead they are:')
    disp(diag(lambda))
end

%--------------------------------------------------------------------------
% EIGENVECTORS
%--------------------------------------------------------------------------
evec0 = double(VecSorted(:,1));         % Assign each eigenvector with the corresponding eigenvalue
evec1 = double(VecSorted(:,2));         % NOTE: in the end, we have not explicitly used these eigenvectors
evec2 = double(VecSorted(:,3));
evec3 = double(VecSorted(:,4));
evec4 = double(VecSorted(:,5));

disp('Time to calculate eigenvalue problem:')
toc
%%
%--------------------------------------------------------------------------
% CONDITIONAL PROBABILITIES
%--------------------------------------------------------------------------
% Evaluate conditional probabilties by substituting in values from above
% and using vpa() to force the simplest form of the output.
tic
disp('...Calculating conditional probabilities')

% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);

% Solve the equations and give an output
% Condition 0: P(t=0) = [P0(0) = 1, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 0]
Psol_0 = dsolve(eqn,[P0(0)==1 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==0]);
P04(t) = real(vpa(Psol_0.P4));
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

%--------------------------------------------------------------------------
% DAMAGE CONTROL
%--------------------------------------------------------------------------
if eval0 == 0
    disp('Conditional probabilities are good to go!')
else
    tic
    digits(10);
    disp('Fixing conditional probabilities because of non-zero eval0')
    
    eval0_string = num2str(eval0);                  % define a string of the non-zero value of eval0
    eval0_str = eval0_string(1:7);                  % only take the first 7 characters of the string
    
    P00_cell = sym2cell(P00);                       % convert P00 expression to a cell
    P00_str = string(P00_cell);                     % then turn it into a string
    P00_str_new = strrep(P00_str,eval0_str,'0*');   % replace the non-zero eval with zero
    P00(t) = str2sym(P00_str_new);                  % turn the string back into a symbolic expression
    
    P01_cell = sym2cell(P01);                  
    P01_str = string(P01_cell);                 
    P01_str_new = strrep(P01_str,eval0_str,'0*');    
    P01(t) = str2sym(P01_str_new);     
    
    P02_cell = sym2cell(P02);                  
    P02_str = string(P02_cell); 
    P02_str_new = strrep(P02_str,eval0_str,'0*');    
    P02(t) = str2sym(P02_str_new);
    
    P03_cell = sym2cell(P03);                  
    P03_str = string(P03_cell);                 
    P03_str_new = strrep(P03_str,eval0_str,'0*');    
    P03(t) = str2sym(P03_str_new);
    
    P04_cell = sym2cell(P04);                  
    P04_str = string(P04_cell);                 
    P04_str_new = strrep(P04_str,eval0_str,'0*');    
    P04(t) = str2sym(P04_str_new);
    
    P10_cell = sym2cell(P10);                  
    P10_str = string(P10_cell);                 
    P10_str_new = strrep(P10_str,eval0_str,'0*');    
    P10(t) = str2sym(P10_str_new);
    
    P11_cell = sym2cell(P11);                  
    P11_str = string(P11_cell);                 
    P11_str_new = strrep(P11_str,eval0_str,'0*');    
    P11(t) = str2sym(P11_str_new);
    
    P12_cell = sym2cell(P12);                  
    P12_str = string(P12_cell);                 
    P12_str_new = strrep(P12_str,eval0_str,'0*');    
    P12(t) = str2sym(P12_str_new);
    
    P13_cell = sym2cell(P13);                  
    P13_str = string(P13_cell);                 
    P13_str_new = strrep(P13_str,eval0_str,'0*');    
    P13(t) = str2sym(P13_str_new);
    
    P14_cell = sym2cell(P14);                  
    P14_str = string(P14_cell);                 
    P14_str_new = strrep(P14_str,eval0_str,'0*');    
    P14(t) = str2sym(P14_str_new);
    
    P20_cell = sym2cell(P20);                  
    P20_str = string(P20_cell);                 
    P20_str_new = strrep(P20_str,eval0_str,'0*');    
    P20(t) = str2sym(P20_str_new);
    
    P21_cell = sym2cell(P21);                  
    P21_str = string(P21_cell);                 
    P21_str_new = strrep(P21_str,eval0_str,'0*');    
    P21(t) = str2sym(P21_str_new);
    
    P22_cell = sym2cell(P22);                  
    P22_str = string(P22_cell);                 
    P22_str_new = strrep(P22_str,eval0_str,'0*');    
    P22(t) = str2sym(P22_str_new);
    
    P23_cell = sym2cell(P23);                  
    P23_str = string(P23_cell);                 
    P23_str_new = strrep(P23_str,eval0_str,'0*');    
    P23(t) = str2sym(P23_str_new);
   
    P24_cell = sym2cell(P24);                  
    P24_str = string(P24_cell);                 
    P24_str_new = strrep(P24_str,eval0_str,'0*');    
    P24(t) = str2sym(P24_str_new);
    
    P30_cell = sym2cell(P30);                  
    P30_str = string(P30_cell);                 
    P30_str_new = strrep(P30_str,eval0_str,'0*');    
    P30(t) = str2sym(P30_str_new);
    
    P31_cell = sym2cell(P31);                  
    P31_str = string(P31_cell);                 
    P31_str_new = strrep(P31_str,eval0_str,'0*');    
    P31(t) = str2sym(P31_str_new);
    
    P32_cell = sym2cell(P32);                  
    P32_str = string(P32_cell);                 
    P32_str_new = strrep(P32_str,eval0_str,'0*');    
    P32(t) = str2sym(P32_str_new);
    
    P33_cell = sym2cell(P33);                  
    P33_str = string(P33_cell);                 
    P33_str_new = strrep(P33_str,eval0_str,'0*');    
    P33(t) = str2sym(P33_str_new);
    
    P34_cell = sym2cell(P34);                  
    P34_str = string(P34_cell);                 
    P34_str_new = strrep(P34_str,eval0_str,'0*');    
    P34(t) = str2sym(P34_str_new);
    
    P40_cell = sym2cell(P40);                  
    P40_str = string(P40_cell);                 
    P40_str_new = strrep(P40_str,eval0_str,'0*');    
    P40(t) = str2sym(P40_str_new);
    
    P41_cell = sym2cell(P41);                  
    P41_str = string(P41_cell);                 
    P41_str_new = strrep(P41_str,eval0_str,'0*');    
    P41(t) = str2sym(P41_str_new);
    
    P42_cell = sym2cell(P42);                  
    P42_str = string(P42_cell);                 
    P42_str_new = strrep(P42_str,eval0_str,'0*');    
    P42(t) = str2sym(P42_str_new);
    
    P43_cell = sym2cell(P43);                  
    P43_str = string(P43_cell);                 
    P43_str_new = strrep(P43_str,eval0_str,'0*');    
    P43(t) = str2sym(P43_str_new);
    
    P44_cell = sym2cell(P44);                  
    P44_str = string(P44_cell);                 
    P44_str_new = strrep(P44_str,eval0_str,'0*');    
    P44(t) = str2sym(P44_str_new);
    
    disp('Fixing conditional probabilities took:')
    toc
end
  

% Define a matrix of the conditional probabilities: (Pi-->j with i is initial condition)
 cP = [P00(t), P10(t), P20(t), P30(t), P40(t);
       P01(t), P11(t), P21(t), P31(t), P41(t);
       P02(t), P12(t), P22(t), P32(t), P42(t);
       P03(t), P13(t), P23(t), P33(t), P43(t);
       P04(t), P14(t), P24(t), P34(t), P44(t)];
% Row = final state & Column = initial condition

disp('Time to calculate the rate equation ODE')
toc

% P00(t) = vpa(subs(P00));    
% P01(t) = vpa(subs(P01));
% P02(t) = vpa(subs(P02));
% P03(t) = vpa(subs(P03));
% 
% P10(t) = vpa(subs(P10));
% P11(t) = vpa(subs(P11));
% P12(t) = vpa(subs(P12));
% P13(t) = vpa(subs(P13));
% 
% P20(t) = vpa(subs(P20));
% P21(t) = vpa(subs(P21));
% P22(t) = vpa(subs(P22));
% P23(t) = vpa(subs(P23));
% 
% P30(t) = vpa(subs(P30));
% P31(t) = vpa(subs(P31));
% P32(t) = vpa(subs(P32));
% P33(t) = vpa(subs(P33));
% Note: these output as symbolic functions ('symfun') need to evaluate at a
% time and call it with double() to turn it into a number.

%%
%--------------------------------------------------------------------------
% EQUILIBRIUM PROBABILITIES
%--------------------------------------------------------------------------
% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P0EQ = double(P00(inf));
P1EQ = double(P11(inf));
P2EQ = double(P22(inf));
P3EQ = double(P33(inf));
P4EQ = double(P44(inf));
 
Peq = [P0EQ; P1EQ; P2EQ; P3EQ; P4EQ];

if sum(Peq) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('PROBLEM: Equilibrium probabilities DO NOT sum to 1.')
    disp('Instead they sum to');
    disp(sum(Peq))
end

%--------------------------------------------------------------------------
% FRET STATES
%--------------------------------------------------------------------------
% Define FRET states in a vector A:
A0 = 0.34;  % These are just guesses for now
A1 = 0.42;
A2 = 0.51;
A3 = 0.67;
A4 = 0.49;

A = [double(A0); double(A1); double(A2); double(A3); double(A4)];
A = double(A);

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


if abs(double(C2(0)) - double(msq)) == 0
    disp('Mean of the square <A^2> matches C2(t=0)!')
elseif abs(double(C2(0)) - double(msq)) <= 1e-7
    disp('Mean of square (msq) <A^2> matches C2(t=0) within 1e-7 - probably rounding')
else 
    disp('PROBLEM: mean of the square <A^2> DOES NOT match C2(t=0)')
    disp('C2(t=0) is')
    disp(C2(0))
    disp('and mean of the square (msq) <A^2> is')
    disp(msq)
end

if double(C2(10^20)) == double(sqm)
    disp('Square of the mean <A>^2 matches C2(t=inf)!')
elseif double(C2(10^20)) < 1e-10 && double(sqm) < 1e-10
    disp('Both C2(t = big) and <A>^2 (sqm) are less than 1e-10, so they match close enough to zero.')
else 
    disp('PROBLEM: square of the mean <A>^2 DOES NOT match C2(t=inf)')
    disp('C2(t=inf) is')
    disp(double(C2(10^20)))
    disp('and <A>^2 (sqm) is')
    disp(double(sqm))
end

%--------------------------------------------------------------------------
% Plot two point TCF
%--------------------------------------------------------------------------

close all
figure
TCF2pt = fplot(C2(t),[0.00001,1]);
title('Analytical Two point TCF','FontSize',18)
xlabel('Time (\tau_1)','FontSize',14);
ylabel('C^{(2)}(\tau)','FontSize',14);
% set(gca,'XScale','log')

 ax = gca;
 ax.XScale = 'log';

%%
%--------------------------------------------------------------------------
% Four point TCF:
%--------------------------------------------------------------------------

% Define times:
Npts = 50;
tau1vec = logspace(0,5,Npts);   % Vector of times
tau2 = 0;                       % Tau2 is a constant
tau3vec = tau1vec;              % Tau3 is the same vector of time as tau1

t1 = tau1vec;
t2 = tau2;
t3 = tau3vec;
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
C4mat = zeros(length(t1));
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

C4diff =double( C4 - (C2(inf))^2 );

disp('Time to calculate the four point TCF (C4):');
toc
%%
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



% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  August 2019
%
% PURPOSE:  Evaluate the seven state (01234) conditional probabilties with a set of 
%           rates then the  2 point TCF and 4 point TCF for that system
%
%       MODEL:  01  02  03  1   1'  2   20
%       CODE:   0   1   2   3   4   5   6
%   
%       NOTE:   We have assumed two loops, between the zero states, and between
%               1', 2 and 20 (so two detailed balance conditions)
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
% MODIFICATIONS: Sept 2019 - calculate C2 for a vector of times
%                            plot C4 for several values of tau2
%
%__________________________________________________________________________

close all
clear all       % CAREFUL: sometimes if the variables are not cleared the eigenvalue 
clc             % calculated from eig() and while solving the ODE will not match.
                % Then all of the subsitutions will break, this is probably
                % what happened if sum(Peq) = 0.
%%
tic
t = sym('t');
syms P0(t) P1(t) P2(t) P3(t) P4(t) P5(t) P6(t)

%--------------------------------------------------------------------------
% Define rates:
%--------------------------------------------------------------------------

rate_str = '1000s10s1s1000'; % This will name the png's of the C2 and C4 figure

k01 = 1000;
k02 = 2000;
            k03 = 0;        % Move rates that are zero out of the way.
            k04 = 0;
            k05 = 0;
            k06 = 0;
k10 = 3000;
k12 = 4000;
k13 = 50;
k14 = 60; 
            k15 = 0;
            k16 = 0;
% k20 = 13; detailed balance
k21 = 7000;
            k23 = 0;
            k24 = 0;
            k25 = 0;
            k26 = 0;
            k30 = 0;
k31 = 80;
            k32 = 0;
            k34 = 0;
            k35 = 0;
            k36 = 0;
            k40 = 0;
k41 = 90;
            k42 = 0;
            k43 = 0;
k45 = 50;
k46 = 8;
            k50 = 0;
            k51 = 0;
            k52 = 0;
            k53 = 0;
k54 = 6;
k56 = 9300;
            k60 = 0;
            k61 = 0;
            k62 = 0;
            k63 = 0;
% k64 = ; detailed balance
k65 = 1400;
% NOTE: I put in the rates we know are zero in by hand. For this model we
% have 10 non-zero rates.


% Detailed balance condition
k20 = (k10 * k02 * k21)/(k01 * k12);
k64 = (k46 * k65 * k54)/(k45 * k56);

%--------------------------------------------------------------------------
% EIGENVALUE PROBLEM
%--------------------------------------------------------------------------
tic
disp('...Solving eigenvalue problem')
% Define rate matrix
K = [-(k01+k02+k03+k04+k05+k06), k10,           k20,              k30,              k40,               k50,               k60;...
        k01,      -(k10+k12+k13+k14+k15+k16),   k21,              k31,              k41,               k51,               k61;...
        k02,             k12,       -(k20+k21+k23+k24+k25+k26),   k32,              k42,               k52,               k62;...
        k03,             k13,                   k23,             -k31,              k43,               k53,               k63;...
        k04,             k14,                   k24,              k34,    -(k40+k41+k42+k43+k45+k46),  k54,               k64;...
        k05,             k15,                   k25,              k35,              k45,      -(k50+k51+k52+k53+k54+k56)  k65;...
        k06,             k16,                   k26,              k36,              k46,               k56,     -(k60+k61+k62+k63+k64+k65)];
   
    
% K(4,4) =  -(k30+k31+k32+k34+k35+36) % I don't know why this sums to -56
% when only k31 is non-zero and it is 20

    % Vector of populations (or probabilities)
P(t) = [P0(t); P1(t); P2(t); P3(t); P4(t);P5(t);P6(t)];


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
eval5 = double(lambdaSorted(6,6));
eval6 = double(lambdaSorted(7,7));

if eval0<=0 && eval1<=0 && eval2 <=0 && eval3 <=0 && eval4 <=0 && eval5<=0 && eval6<=0      % For the exponentials to decay, need all eigenvalues to be < 0.
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
evec5 = double(VecSorted(:,6));
evec6 = double(VecSorted(:,7));

elapsedTime = toc;
task_str = 'to calculate eigenvalue problem';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);

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
% Condition 0: P(t=0) = [P0(0) = 1, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 0, P5(0) = 0, P6(0) = 0]
tic
Psol_0 = dsolve(eqn,[P0(0)==1 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==0 P5(0)==0 P6(0)==0]);
P04(t) = real(vpa(Psol_0.P4));
P00(t) = real(vpa(Psol_0.P0));
P01(t) = real(vpa(Psol_0.P1));
P02(t) = real(vpa(Psol_0.P2));
P03(t) = real(vpa(Psol_0.P3));
P04(t) = real(vpa(Psol_0.P4));
P05(t) = real(vpa(Psol_0.P5));
P06(t) = real(vpa(Psol_0.P6));
toc

% Condition 1: P(t=0) = [P0(0) = 0, P1(0) = 1, P2(0) = 0, P3(0) = 0, P4(0) = 0, P5(0) = 0, P6(0) = 0]
tic
Psol_1 = dsolve(eqn,[P0(0) == 0 P1(0)==1 P2(0)==0 P3(0)==0 P4(0)==0 P5(0)==0 P6(0)==0]);
P10(t) = real(vpa(Psol_1.P0));
P11(t) = real(vpa(Psol_1.P1));
P12(t) = real(vpa(Psol_1.P2));
P13(t) = real(vpa(Psol_1.P3));
P14(t) = real(vpa(Psol_1.P4));
P15(t) = real(vpa(Psol_1.P5));
P16(t) = real(vpa(Psol_1.P6));
toc

% Condition 2: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 1, P3(0) = 0, P4(0) = 0, P5(0) = 0, P6(0) = 0]
tic
Psol_2 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==1 P3(0)==0 P4(0)==0 P5(0)==0 P6(0)==0]);
P20(t) = real(vpa(Psol_2.P0));
P21(t) = real(vpa(Psol_2.P1));
P22(t) = real(vpa(Psol_2.P2));
P23(t) = real(vpa(Psol_2.P3));
P24(t) = real(vpa(Psol_2.P4));
P25(t) = real(vpa(Psol_2.P5));
P26(t) = real(vpa(Psol_2.P6));
toc

% Condition 3: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 1, P4(0) = 0, P5(0) = 0, P6(0) = 0]
tic
Psol_3 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==1 P4(0)==0 P5(0)==0 P6(0)==0]);
P30(t) = real(vpa(Psol_3.P0));
P31(t) = real(vpa(Psol_3.P1));
P32(t) = real(vpa(Psol_3.P2));
P33(t) = real(vpa(Psol_3.P3));
P34(t) = real(vpa(Psol_3.P4));
P35(t) = real(vpa(Psol_3.P5));
P36(t) = real(vpa(Psol_3.P6));
toc

% Condition 4: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 1, P5(0) = 0, P6(0) = 0]
Psol_4 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==1 P5(0)==0 P6(0)==0]);
P40(t) = real(vpa(Psol_4.P0));
P41(t) = real(vpa(Psol_4.P1));
P42(t) = real(vpa(Psol_4.P2));
P43(t) = real(vpa(Psol_4.P3));
P44(t) = real(vpa(Psol_4.P4));
P45(t) = real(vpa(Psol_4.P5));
P46(t) = real(vpa(Psol_4.P6));

% Condition 5: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 0, P5(0) = 1, P6(0) = 0]
tic
Psol_5 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==0 P5(0)==1 P6(0)==0]);
P50(t) = real(vpa(Psol_5.P0));
P51(t) = real(vpa(Psol_5.P1));
P52(t) = real(vpa(Psol_5.P2));
P53(t) = real(vpa(Psol_5.P3));
P54(t) = real(vpa(Psol_5.P4));
P55(t) = real(vpa(Psol_5.P5));
P56(t) = real(vpa(Psol_5.P6));
toc

% Condition 6: P(t=0) = [P0(0) = 0, P1(0) = 0, P2(0) = 0, P3(0) = 0, P4(0) = 0, P5(0) = 0, P6(0) = 1]
tic
Psol_6 = dsolve(eqn,[P0(0) == 0 P1(0)==0 P2(0)==0 P3(0)==0 P4(0)==0 P5(0)==0 P6(0)==1]);
P60(t) = real(vpa(Psol_6.P0));
P61(t) = real(vpa(Psol_6.P1));
P62(t) = real(vpa(Psol_6.P2));
P63(t) = real(vpa(Psol_6.P3));
P64(t) = real(vpa(Psol_6.P4));
P65(t) = real(vpa(Psol_6.P5));
P66(t) = real(vpa(Psol_6.P6));
toc

disp('Time to calculate conditional probabilities')
toc

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
    eval0_str = eval0_string(1:6);                  % only take the first 7 characters of the string
    
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
    
    P05_cell = sym2cell(P05);                  
    P05_str = string(P05_cell);                 
    P05_str_new = strrep(P05_str,eval0_str,'0*');    
    P05(t) = str2sym(P05_str_new);
    
    P06_cell = sym2cell(P06);                  
    P06_str = string(P06_cell);                 
    P06_str_new = strrep(P06_str,eval0_str,'0*');    
    P06(t) = str2sym(P06_str_new);
    
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
    
    P15_cell = sym2cell(P15);                  
    P15_str = string(P15_cell);                 
    P15_str_new = strrep(P15_str,eval0_str,'0*');    
    P15(t) = str2sym(P15_str_new);
    
    P16_cell = sym2cell(P16);                  
    P16_str = string(P16_cell);                 
    P16_str_new = strrep(P16_str,eval0_str,'0*');    
    P16(t) = str2sym(P16_str_new);
    
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
    
    P25_cell = sym2cell(P25);                  
    P25_str = string(P25_cell);                 
    P25_str_new = strrep(P25_str,eval0_str,'0*');    
    P25(t) = str2sym(P25_str_new);
   
    P26_cell = sym2cell(P26);                  
    P26_str = string(P26_cell);                 
    P26_str_new = strrep(P26_str,eval0_str,'0*');    
    P26(t) = str2sym(P26_str_new);
    
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
    
    P35_cell = sym2cell(P35);                  
    P35_str = string(P35_cell);                 
    P35_str_new = strrep(P35_str,eval0_str,'0*');    
    P35(t) = str2sym(P35_str_new);
    
    P36_cell = sym2cell(P36);                  
    P36_str = string(P36_cell);                 
    P36_str_new = strrep(P36_str,eval0_str,'0*');    
    P36(t) = str2sym(P36_str_new);
    
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
    
    P45_cell = sym2cell(P45);                  
    P45_str = string(P45_cell);                 
    P45_str_new = strrep(P45_str,eval0_str,'0*');    
    P45(t) = str2sym(P45_str_new);
    
    P46_cell = sym2cell(P46);                  
    P46_str = string(P46_cell);                 
    P46_str_new = strrep(P46_str,eval0_str,'0*');    
    P46(t) = str2sym(P46_str_new);
    
    P50_cell = sym2cell(P50);                  
    P50_str = string(P50_cell);                 
    P50_str_new = strrep(P50_str,eval0_str,'0*');    
    P50(t) = str2sym(P50_str_new);
    
    P51_cell = sym2cell(P51);                  
    P51_str = string(P51_cell);                 
    P51_str_new = strrep(P51_str,eval0_str,'0*');    
    P51(t) = str2sym(P51_str_new);
    
    P52_cell = sym2cell(P52);                  
    P52_str = string(P52_cell);                 
    P52_str_new = strrep(P52_str,eval0_str,'0*');    
    P52(t) = str2sym(P52_str_new);
    
    P53_cell = sym2cell(P53);                  
    P53_str = string(P53_cell);                 
    P53_str_new = strrep(P53_str,eval0_str,'0*');    
    P53(t) = str2sym(P53_str_new);
    
    P54_cell = sym2cell(P54);                  
    P54_str = string(P54_cell);                 
    P54_str_new = strrep(P54_str,eval0_str,'0*');    
    P54(t) = str2sym(P54_str_new);
    
    P55_cell = sym2cell(P55);                  
    P55_str = string(P55_cell);                 
    P55_str_new = strrep(P55_str,eval0_str,'0*');    
    P55(t) = str2sym(P55_str_new);
    
    P56_cell = sym2cell(P56);                  
    P56_str = string(P56_cell);                 
    P56_str_new = strrep(P56_str,eval0_str,'0*');    
    P56(t) = str2sym(P56_str_new);
    
    P60_cell = sym2cell(P60);                  
    P60_str = string(P60_cell);                 
    P60_str_new = strrep(P60_str,eval0_str,'0*');    
    P60(t) = str2sym(P60_str_new);
    
    P61_cell = sym2cell(P61);                  
    P61_str = string(P61_cell);                 
    P61_str_new = strrep(P61_str,eval0_str,'0*');    
    P61(t) = str2sym(P61_str_new);
    
    P62_cell = sym2cell(P62);                  
    P62_str = string(P62_cell);                 
    P62_str_new = strrep(P62_str,eval0_str,'0*');    
    P62(t) = str2sym(P62_str_new);
    
    P63_cell = sym2cell(P63);                  
    P63_str = string(P63_cell);                 
    P63_str_new = strrep(P63_str,eval0_str,'0*');    
    P63(t) = str2sym(P63_str_new);
    
    P64_cell = sym2cell(P64);                  
    P64_str = string(P64_cell);                 
    P64_str_new = strrep(P64_str,eval0_str,'0*');    
    P64(t) = str2sym(P64_str_new);
    
    P65_cell = sym2cell(P65);                  
    P65_str = string(P65_cell);                 
    P65_str_new = strrep(P65_str,eval0_str,'0*');    
    P65(t) = str2sym(P65_str_new);
    
    P66_cell = sym2cell(P66);                  
    P66_str = string(P66_cell);                 
    P66_str_new = strrep(P66_str,eval0_str,'0*');    
    P66(t) = str2sym(P66_str_new);
    
elapsedTime = toc;
task_str = 'to fix the conditional probabilities:';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);

end
  

% Define a matrix of the conditional probabilities: (Pi-->j with i is initial condition)
 cP = [P00(t), P10(t), P20(t), P30(t), P40(t), P50(t), P60(t);
       P01(t), P11(t), P21(t), P31(t), P41(t), P51(t), P61(t);
       P02(t), P12(t), P22(t), P32(t), P42(t), P52(t), P62(t);
       P03(t), P13(t), P23(t), P33(t), P43(t), P53(t), P63(t);
       P04(t), P14(t), P24(t), P34(t), P44(t), P54(t), P64(t);
       P05(t), P15(t), P25(t), P35(t), P45(t), P55(t), P65(t);
       P06(t), P16(t), P26(t), P36(t), P46(t), P56(t), P66(t)];
% Row = final state & Column = initial condition

disp('Time to calculate the rate equation ODE')
toc


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
P5EQ = double(P55(inf));
P6EQ = double(P66(inf));
 
Peq = [P0EQ; P1EQ; P2EQ; P3EQ; P4EQ; P5EQ; P6EQ];

if double(sum(Peq)) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('Equilibrium probabilities sum to:')
    disp(sum(Peq))
end

%--------------------------------------------------------------------------
% FRET STATES
%--------------------------------------------------------------------------
% Define FRET states in a vector A:
A0 = 0.75;   % These are just guesses for now
A1 = 0.75;
A2 = 0.75;   % We are saying the 0 states, 1 states, and 2 states are degenerate.
A3 = 0.55;
A4 = 0.55;
A5 = 0.3;
A6 = 0.3;

A = [double(A0); double(A1); double(A2); double(A3); double(A4); double(A5); double(A6)];
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

C2sym(t) = 0*t;
for i = 1:numel(A)
    for j = 1:numel(A)
        
        C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
        C2sym(t) = C2sym(t) + C2temp(t);
    end
end


if abs(double(C2sym(0)) - double(msq)) == 0
    disp('Mean of the square <A^2> matches C2(t=0)!')
elseif abs(double(C2sym(0)) - double(msq)) <= 1e-7
    disp('Mean of square (msq) <A^2> matches C2(t=0) within 1e-7 - probably rounding')
else 
    disp('PROBLEM: mean of the square <A^2> DOES NOT match C2(t=0)')
    disp('C2(t=0) is')
    disp(C2sym(0))
    disp('and mean of the square (msq) <A^2> is')
    disp(msq)
end

if double(C2sym(10^20)) == double(sqm)
    disp('Square of the mean <A>^2 matches C2(t=inf)!')
elseif double(C2sym(10^20)) < 1e-10 && double(sqm) < 1e-10
    disp('Both C2(t = big) and <A>^2 (sqm) are less than 1e-10, so they match close enough to zero.')
else 
    disp('PROBLEM: square of the mean <A>^2 DOES NOT match C2(t=inf)')
    disp('C2(t=inf) is')
    disp(double(C2sym(10^20)))
    disp('and <A>^2 (sqm) is')
    disp(double(sqm))
end


%--------------------------------------------------------------------------
% Evaluate C2 over a range of t's
%--------------------------------------------------------------------------

Npts = 150;
timeArray = [1:9,logspace(1,6.4771212,Npts)];
C2 = double(C2sym(timeArray));



%--------------------------------------------------------------------------
% Plot two point TCF
%--------------------------------------------------------------------------

%close all
figure

%subplot(1,2,1)
%TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
TCF2pt = plot(C2,'LineWidth',2);

title('Analytical Two point TCF','FontSize',18)
xlabel('Time (\tau_1)','FontSize',14);
ylabel('C^{(2)}(\tau)','FontSize',14);
%xlim([10^-5 10^0])
% ylim([C2(500) C2(0)])

ax = gca;
ax.XScale = 'log';
set(gca,'yscale','log')

saveName = ['C2_',rate_str];
saveas(TCF2pt,saveName, 'png')

%%
%--------------------------------------------------------------------------
% Four point TCF:
%--------------------------------------------------------------------------

% Define times:
Npts = 100;
tau1vec = logspace(0,5,Npts);   % Vector of times
%tau1vec = timeArray;
tau2 = [0.0001,0.1, 0.5, 1, 2];                       % Tau2 is a constant
tau3vec = tau1vec;              % Tau3 is the same vector of time as tau1

t1 = tau1vec;
t2 = tau2;
t3 = tau3vec;
% These are what get inserted into the conditional probability functions.


tic
disp('... Calculating the Conditional Probability Matrix');
% Initialize size of 3D array for conditional probability for t1:
cP_t1 = zeros(5,5,length(t1));

% Define the vector in each position
cP_t1(1,1,:) = double(P00(t1));
cP_t1(2,1,:) = double(P01(t1));
cP_t1(3,1,:) = double(P02(t1));
cP_t1(4,1,:) = double(P03(t1));
cP_t1(5,1,:) = double(P04(t1));
cP_t1(6,1,:) = double(P05(t1));
cP_t1(7,1,:) = double(P06(t1));

cP_t1(1,2,:) = double(P10(t1));
cP_t1(2,2,:) = double(P11(t1));
cP_t1(3,2,:) = double(P12(t1));
cP_t1(4,2,:) = double(P13(t1));
cP_t1(5,2,:) = double(P14(t1));
cP_t1(6,2,:) = double(P15(t1));
cP_t1(7,2,:) = double(P16(t1));

cP_t1(1,3,:) = double(P20(t1));
cP_t1(2,3,:) = double(P21(t1));
cP_t1(3,3,:) = double(P22(t1));
cP_t1(4,3,:) = double(P23(t1));
cP_t1(5,3,:) = double(P24(t1));
cP_t1(6,3,:) = double(P25(t1));
cP_t1(7,3,:) = double(P26(t1));

cP_t1(1,4,:) = double(P30(t1));
cP_t1(2,4,:) = double(P31(t1));
cP_t1(3,4,:) = double(P32(t1));
cP_t1(4,4,:) = double(P33(t1));
cP_t1(5,4,:) = double(P34(t1));
cP_t1(6,4,:) = double(P35(t1));
cP_t1(7,4,:) = double(P36(t1));

cP_t1(1,5,:) = double(P40(t1));
cP_t1(2,5,:) = double(P41(t1));
cP_t1(3,5,:) = double(P42(t1));
cP_t1(4,5,:) = double(P43(t1));
cP_t1(5,5,:) = double(P44(t1));
cP_t1(6,5,:) = double(P45(t1));
cP_t1(7,5,:) = double(P46(t1));

cP_t1(1,6,:) = double(P50(t1));
cP_t1(2,6,:) = double(P51(t1));
cP_t1(3,6,:) = double(P52(t1));
cP_t1(4,6,:) = double(P53(t1));
cP_t1(5,6,:) = double(P54(t1));
cP_t1(6,6,:) = double(P55(t1));
cP_t1(7,6,:) = double(P56(t1));

cP_t1(1,7,:) = double(P60(t1));
cP_t1(2,7,:) = double(P61(t1));
cP_t1(3,7,:) = double(P62(t1));
cP_t1(4,7,:) = double(P63(t1));
cP_t1(5,7,:) = double(P64(t1));
cP_t1(6,7,:) = double(P65(t1));
cP_t1(7,7,:) = double(P66(t1));

% Initialize size of 3D array for conditional probability for t2:
cP_t2 = zeros(4,4,length(t2));
% Define the vector in each position
cP_t2(1,1,:) = double(P00(t2));
cP_t2(2,1,:) = double(P01(t2));
cP_t2(3,1,:) = double(P02(t2));
cP_t2(4,1,:) = double(P03(t2));
cP_t2(5,1,:) = double(P04(t2));
cP_t2(6,1,:) = double(P05(t2));
cP_t2(7,1,:) = double(P06(t2));

cP_t2(1,2,:) = double(P10(t2));
cP_t2(2,2,:) = double(P11(t2));
cP_t2(3,2,:) = double(P12(t2));
cP_t2(4,2,:) = double(P13(t2));
cP_t2(5,2,:) = double(P14(t2));
cP_t2(6,2,:) = double(P15(t2));
cP_t2(7,2,:) = double(P16(t2));

cP_t2(1,3,:) = double(P20(t2));
cP_t2(2,3,:) = double(P21(t2));
cP_t2(3,3,:) = double(P22(t2));
cP_t2(4,3,:) = double(P23(t2));
cP_t2(5,3,:) = double(P24(t2));
cP_t2(6,3,:) = double(P25(t2));
cP_t2(7,3,:) = double(P26(t2));

cP_t2(1,4,:) = double(P30(t2));
cP_t2(2,4,:) = double(P31(t2));
cP_t2(3,4,:) = double(P32(t2));
cP_t2(4,4,:) = double(P33(t2));
cP_t2(5,4,:) = double(P34(t2));
cP_t2(6,4,:) = double(P35(t2));
cP_t2(7,4,:) = double(P36(t2));

cP_t2(1,5,:) = double(P40(t2));
cP_t2(2,5,:) = double(P41(t2));
cP_t2(3,5,:) = double(P42(t2));
cP_t2(4,5,:) = double(P43(t2));
cP_t2(5,5,:) = double(P44(t2));
cP_t2(6,5,:) = double(P45(t2));
cP_t2(7,5,:) = double(P46(t2));

cP_t2(1,6,:) = double(P50(t2));
cP_t2(2,6,:) = double(P51(t2));
cP_t2(3,6,:) = double(P52(t2));
cP_t2(4,6,:) = double(P53(t2));
cP_t2(5,6,:) = double(P54(t2));
cP_t2(6,6,:) = double(P55(t2));
cP_t2(7,6,:) = double(P56(t2));

cP_t2(1,7,:) = double(P60(t2));
cP_t2(2,7,:) = double(P61(t2));
cP_t2(3,7,:) = double(P62(t2));
cP_t2(4,7,:) = double(P63(t2));
cP_t2(5,7,:) = double(P64(t2));
cP_t2(6,7,:) = double(P65(t2));
cP_t2(7,7,:) = double(P66(t2));

% Initialize size of 3D array for conditional probability for t3:
cP_t3 = zeros(5,5,length(t3));
% Define the vector in each position
cP_t3(1,1,:) = double(P00(t3));
cP_t3(2,1,:) = double(P01(t3));
cP_t3(3,1,:) = double(P02(t3));
cP_t3(4,1,:) = double(P03(t3));
cP_t3(5,1,:) = double(P04(t3));
cP_t3(6,1,:) = double(P05(t3));
cP_t3(7,1,:) = double(P06(t3));

cP_t3(1,2,:) = double(P10(t3));
cP_t3(2,2,:) = double(P11(t3));
cP_t3(3,2,:) = double(P12(t3));
cP_t3(4,2,:) = double(P13(t3));
cP_t3(5,2,:) = double(P14(t3));
cP_t3(6,2,:) = double(P15(t3));
cP_t3(7,2,:) = double(P16(t3));

cP_t3(1,3,:) = double(P20(t3));
cP_t3(2,3,:) = double(P21(t3));
cP_t3(3,3,:) = double(P22(t3));
cP_t3(4,3,:) = double(P23(t3));
cP_t3(5,3,:) = double(P24(t3));
cP_t3(6,3,:) = double(P25(t3));
cP_t3(7,3,:) = double(P26(t3));

cP_t3(1,4,:) = double(P30(t3));
cP_t3(2,4,:) = double(P31(t3));
cP_t3(3,4,:) = double(P32(t3));
cP_t3(4,4,:) = double(P33(t3));
cP_t3(5,4,:) = double(P34(t3));
cP_t3(6,4,:) = double(P35(t3));
cP_t3(7,4,:) = double(P36(t3));

cP_t3(1,5,:) = double(P40(t3));
cP_t3(2,5,:) = double(P41(t3));
cP_t3(3,5,:) = double(P42(t3));
cP_t3(4,5,:) = double(P43(t3));
cP_t3(5,5,:) = double(P44(t3));
cP_t3(6,5,:) = double(P45(t3));
cP_t3(7,5,:) = double(P46(t3));

cP_t3(1,6,:) = double(P50(t3));
cP_t3(2,6,:) = double(P51(t3));
cP_t3(3,6,:) = double(P52(t3));
cP_t3(4,6,:) = double(P53(t3));
cP_t3(5,6,:) = double(P54(t3));
cP_t3(6,6,:) = double(P55(t3));
cP_t3(7,6,:) = double(P56(t3));

cP_t3(1,7,:) = double(P60(t3));
cP_t3(2,7,:) = double(P61(t3));
cP_t3(3,7,:) = double(P62(t3));
cP_t3(4,7,:) = double(P63(t3));
cP_t3(5,7,:) = double(P64(t3));
cP_t3(6,7,:) = double(P65(t3));
cP_t3(7,7,:) = double(P66(t3));

elapsedTime = toc;
task_str = ' to calculate Conditional Probability Matrix';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);

%%

disp('... Calculating the 4 point TCF');
tic
C4mat = zeros(length(t1));
 %-------------------------------------------------------------------------
 % Iterate over all the Combinations of FRET States
 %-------------------------------------------------------------------------
 figure         % Create figure OUTSIDE loop so plots overlay
 
 for g = 1:numel(t2)    % Calculate for values of tau2 = [0.0001,0.1, 0.5, 1, 2]
     
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

C4diff = double( C4 - (C2sym(inf))^2 );

elapsedTime = toc;
task_str = ' to calculate the four point TCF (C4)';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);


%--------------------------------------------------------------------------
% Plot surface of C4
%--------------------------------------------------------------------------
%subplot(1,2,2)
TCF4pt = surf(tau1vec, tau3vec, C4);
hold on
title('Analytical Four-point TCF: C^{(4)}','FontSize',18)
xlabel('Time (\tau_1)','FontSize',14);
ylabel('Time (\tau_3)','FontSize',14);
zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
xlim([10^0 10^1])
ylim([10^0 10^1])
%shading interp

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

 end

view(45,45)
 hold off
 
 colorbar
 % save png of figure
saveName = ['C4_',rate_str];
saveas(TCF4pt,saveName, 'png')
 
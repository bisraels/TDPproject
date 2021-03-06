% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: August 2019
%
% PURPOSE:  Evaluate the four state conditional probabilties with a set of 
%           rates then the  2 point TCF and 4 point TCF for that system
%
% INPUT:    conditional probabilities from FourStateODE_0123 
%           loaded from: 'symCondProb_4state0123_fixed.mat'
%           
%
% MODIFICATIONS:  August 2019 - Added code to calculate two point and four point TCFs
%
% TO FIX: (1) Time step needs to be replaced with the value
%
%--------------------------------------------------------------------------

%Output created by FourStateODE_0123_fixed.m
disp('Loading the conditional Probabilities as a function of rates');
tic
%load('symCondProb_4state0123.mat')        % This is from the old code
load('symCondProb_4state0123_fixed.mat')   % This is the output of the fixed code.
toc

disp('
tic
% syms t
t = sym('t');

k01 = 0.001;
%k02 = 2;
k03 = 0;
k10 = 0.000034;
k12 = 0.44;
k13 = 0;
k20 = 0.153;
k21 = 0.0060;
k23 = 700;
k30 = 0;
k31 = 0;
k32 = 87;
% NOTE: I put in the rates we know are zero in by hand.

% Detailed balance condition
% k20 = (k10 * k02 * k21)/(k10 * k12);
k02 = (k01 * k12 * k20)/(k10 * k21);

% Evaluate the eigenvalues in terms of the rates defined above - produce as
% doubles
eval0 = double(vpa(subs(eval0)));
eval1 = double(vpa(subs(eval1)));
eval2 = double(vpa(subs(eval2)));
eval3 = double(vpa(subs(eval3)));

% Evaluate conditional probabilties by substituting in values from above
% and using vpa() to force the simplest form of the output.

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

% NOTE: When using these expressions, when you plug in a value of t, the
% expression will still be a of class 'sym'.
% To evaluate these expressions as doubles use the following:
%       If t = 0:
%       P00_t1 = double(P00(0))
% Now P00_t1 will be a double


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CALCULATE TCF's     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need:
%       - Pji(t) calcualted above
%       - Pi_eq equilibruim populations (from TDP_probcalc)
%           -> do not need TDP_probcalc - can calc Peq from condProb's
%       - Ai values of FRET states

%--------------------------------------------------------------------------
% Two point TCF:
%--------------------------------------------------------------------------
% syms A0 A1 A2 A3
% syms t        % This forces a static workspace
t = sym('t');   % This allows for a dynamic workspace

% This is a vector for the FRET values - assign values
A0 = 0.34;  % These are just guesses for now
A1 = 0.42;
A2 = 0.51;
A3 = 0.67;

A = [A0; A1; A2; A3];

% Matrix of conditional probabilities Pi-->j with i is initial condition
% cP = [P00(t), P01(t), P02(t), P03(t);...
%      P10(t), P11(t), P12(t), P13(t);...
%      P20(t), P21(t), P22(t), P23(t);...
%      P30(t), P31(t), P32(t), P33(t)];
 cP = [P00(t), P10(t), P20(t), P30(t);...
       P01(t), P11(t), P21(t), P31(t);...
       P02(t), P12(t), P22(t), P32(t);...
       P03(t), P13(t), P23(t), P33(t)];
% Row = final state & Column = initial condition

%%
% load('sym_prob.mat')
% This gives the equilibrium probabilties
% for this state we are using:
% Pi_Cn4 where i = 0, 1, 2, 3
 
% P0_eq = vpa(subs(P0_Cn4));
% P1_eq = vpa(subs(P1_Cn4));
% P2_eq = vpa(subs(P2_Cn4));
% P3_eq = vpa(subs(P3_Cn4));
% 
% Peq = [P0_eq; P1_eq; P2_eq; P3_eq];

% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P0EQ = P00(inf);
P1EQ = P11(inf);
P2EQ = P22(inf);
P3EQ = P33(inf);
 
Peq = [P0EQ; P1EQ; P2EQ; P3EQ];

if sum(Peq) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('Problem: Equilibrium probabilities DO NOT sum to 1.')
end

%%
Amean = sum(A.*Peq);
A = A - Amean;

msq = sum((A.^2).*Peq);     % square of mean <A^2>
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2


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
% ax = gca;
% ax.XScale = 'log';



%%
%--------------------------------------------------------------------------
% Four point TCF:
%--------------------------------------------------------------------------

 

% syms A0 A1 A2 A3
% sym t
t = sym('t');

% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P0EQ = double(P00(inf));
P1EQ = double(P11(inf));
P2EQ = double(P22(inf));
P3EQ = double(P33(inf));

% global Peq
Peq = [P0EQ; P1EQ; P2EQ; P3EQ];

if sum(Peq) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('Problem: Equilibrium probabilities DO NOT sum to 1.')
end

% This is a vector for the FRET values - assign values
A0 = 0.34;  % These are just guesses for now
A1 = 0.42;
A2 = 0.51;
A3 = 0.67;

% global A
A = [double(A0); double(A1); double(A2); double(A3)];

Amean = sum(A.*Peq);
A = A - Amean;

msq = sum((A.^2).*Peq);     % square of mean <A^2>
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2


% syms tau1 tau2 tau3
% tau1 = sym('tau1');
% tau2 = sym('tau2');
% tau3 = sym('tau3');

% t0 = 0;

% Set tau 2
tau2 = 1;

% Create a vector of tau1 and tau3 values
Npts = 50;
tau1vec = logspace(0,3,Npts); % logspace(starting exponent, final exponent, number of bins)
tau3vec = logspace(0,3,Npts);
%tau1vec = [1, 2, 3, 4]; 
%tau3vec = [1, 2, 3, 4];


% Create a matrix to hold the C4's calculated for each (tau1,tau3) pair for
% a set tau 2
C4mat = zeros(length(tau1vec),length(tau3vec));

% Loop over values for tau 1 and tau 3
% for  g = 1:length(tau1vec)
%     tau1 = tau1vec(g);
% for h = 1:length(tau3vec)
%         tau3 = tau3vec(h);

% Or define tau1 and tau3 symbolically and plot the functions
% tau1 = sym('tau1');
% tau2 = sym('tau2');

t1 = tau1vec;
t2 = tau2;
t3 = tau3vec';
        

tic
disp('... Calculating the Conditional Probabilities');
% Initialize size of 3D array for conditional probability for t1:
cP_t1 = zeros(4,4,length(t1));

% Define the vector in each position
cP_t1(1,1,:) = double(P00(t1));
cP_t1(2,1,:) = double(P01(t1));
cP_t1(3,1,:) = double(P02(t1));
cP_t1(4,1,:) = double(P03(t1));

cP_t1(1,2,:) = double(P10(t1));
cP_t1(2,2,:) = double(P11(t1));
cP_t1(3,2,:) = double(P12(t1));
cP_t1(4,2,:) = double(P13(t1));

cP_t1(1,3,:) = double(P20(t1));
cP_t1(2,3,:) = double(P21(t1));
cP_t1(3,3,:) = double(P22(t1));
cP_t1(4,3,:) = double(P23(t1));

cP_t1(1,4,:) = double(P30(t1));
cP_t1(2,4,:) = double(P31(t1));
cP_t1(3,4,:) = double(P32(t1));
cP_t1(4,4,:) = double(P33(t1));



% Initialize size of 3D array for conditional probability for t2:
cP_t2 = zeros(4,4,length(t2));
% Define the vector in each position
cP_t2(1,1,:) = double(P00(t2));
cP_t2(2,1,:) = double(P01(t2));
cP_t2(3,1,:) = double(P02(t2));
cP_t2(4,1,:) = double(P03(t2));

cP_t2(1,2,:) = double(P10(t2));
cP_t2(2,2,:) = double(P11(t2));
cP_t2(3,2,:) = double(P12(t2));
cP_t2(4,2,:) = double(P13(t2));

cP_t2(1,3,:) = double(P20(t2));
cP_t2(2,3,:) = double(P21(t2));
cP_t2(3,3,:) = double(P22(t2));
cP_t2(4,3,:) = double(P23(t2));

cP_t2(1,4,:) = double(P30(t2));
cP_t2(2,4,:) = double(P31(t2));
cP_t2(3,4,:) = double(P32(t2));
cP_t2(4,4,:) = double(P33(t2));

         
% Initialize size of 3D array for conditional probability for t3:
cP_t3 = zeros(4,4,length(t3));
% Define the vector in each position
cP_t3(1,1,:) = double(P00(t3));
cP_t3(2,1,:) = double(P01(t3));
cP_t3(3,1,:) = double(P02(t3));
cP_t3(4,1,:) = double(P03(t3));

cP_t3(1,2,:) = double(P10(t3));
cP_t3(2,2,:) = double(P11(t3));
cP_t3(3,2,:) = double(P12(t3));
cP_t3(4,2,:) = double(P13(t3));

cP_t3(1,3,:) = double(P20(t3));
cP_t3(2,3,:) = double(P21(t3));
cP_t3(3,3,:) = double(P22(t3));
cP_t3(4,3,:) = double(P23(t3));

cP_t3(1,4,:) = double(P30(t3));
cP_t3(2,4,:) = double(P31(t3));
cP_t3(3,4,:) = double(P32(t3));
cP_t3(4,4,:) = double(P33(t3));  
       
disp('Time to calculate Conditional Probabilities:');
      toc 
      
 %C4 = @(tau1, tau2, tau3) 0 * tau1 * tau2 * tau3;
 C4vec = [];
 
 % C4term_val = [];
 % C4mat = zeros(length(A),length(A),length(t1),length(t2));

disp('... Calculating the 4 point TCF');
tic
%-------------------------------------------------------------------------
 % Loop1: Iterate over the time dimension
 %-------------------------------------------------------------------------
%  for timeStep = 1:length(t1)
% 
%      for i = 1:numel(A)
%          for j = 1:numel(A)
%              for k = 1:numel(A)
%                  for l = 1:numel(A)
%                      
%                      % C4temp = @(tau1, tau2, tau3) A(l) * cP_t3(l,k) * A(k) * cP_t2(k,j) * A(j) * cP_t1(j,i) * A(i) * Peq(i);
%                      %C4 = @(tau1, tau2, tau3) C4 + C4temp;
%                      % C4term_val = C4term(tau1, tau2, tau3, i, j, k, l);
% % C4term_val = C4term(t1, t2, t3, A, cP_t1, cP_t2, cP_t3, Peq, i, j, k, l, timeStep);
%                      C4term_val =  A(l) .* cP_t3(l,k,:) .* A(k) .* cP_t2(k,j) .* A(j) .* cP_t1(j,i,timeStep) .* A(i) .* Peq(i);
%                      % C4vec = vertcat(C4vec, C4term_val);
%                      
%                      %C4vec_temp(:,m) = vertcat(C4term_val, C4term_val);
% 
%                      C4vec = [C4vec; C4term_val];
%                      
%                  end
%              end
%          end
%      end
%     
%  end

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


% Plot surface of C4
figure
TCF4point = surf(tau1vec, tau3vec, C4);
title('Analytical Four-point TCF: C^{(4)}','FontSize',18)
xlabel('Time (\tau_1)','FontSize',14);
ylabel('Time (\tau_3)','FontSize',14);
zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';


function C4term_val = C4term(t1, t2, t3, A, cP_t1, cP_t2, cP_t3, Peq, i, j, k, l, timeStep)
C4term_val =  A(l) .* cP_t3(l,k,:) .* A(k) .* cP_t2(k,j) .* A(j) .* cP_t1(j,i,timeStep) .* A(i) .* Peq(i);
end



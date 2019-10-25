% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: August 2019
%
% PURPOSE:  Evaluate the Linear 3-State (123) conditional probabilties with a set of
%           rates then the  2 point TCF and 4 point TCF for that system
%
% INPUT: (1) conditional probabilities from ODE solver: symCondProb_3state123_linear.mat
%
% OUTPUT: (1) Two-point TCF C2
%         (2) Four-point TCF C4
%
% MODIFICATIONS: Adapted from cP2TCF_0123.m
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
verbose_mode = 1; %Set to 1 to see alot of progress updates and print off.


%Output created by ODE solver
disp('Loading the conditional Probabilities as a function of rates');
tic
load('symCondProb_3state123_linear.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')

%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'load the conditional probabilities';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);


tic
% disp('Calculating conditional probabilities using the rates defined')
t = sym('t');

k12 = 1/100;
k21 = 1/300;
k23 = 1/50;
k32 = 1/200;
k13 = 0;
k31 = 0;

% Detailed balance condition
%No detailed balnce condition for the 3-state linear model

% Evaluate the eigenvalues in terms of the rates defined above - produce as doubles
%subs(s) returns a copy of s, replacing symbolic variables in s, with their
%values obtained from the calling function and the MATLAB® Workspace,
% and then evaluates s. Variables with no assigned values remain as variables.
eval1 = double(vpa(subs(eval1)));
eval2 = double(vpa(subs(eval2)));
eval3 = double(vpa(subs(eval3)));

% Evaluate conditional probabilties by substituting in values from above
% and using vpa() to force the simplest form of the output.

%Pij(t) is prob from i--> j, assuming you start in state i: Pi(t=0)=100%=1

P11(t) = vpa(subs(P11));
P12(t) = vpa(subs(P12));
P13(t) = vpa(subs(P13));

P21(t) = vpa(subs(P21));
P22(t) = vpa(subs(P22));
P23(t) = vpa(subs(P23));

P31(t) = vpa(subs(P31));
P32(t) = vpa(subs(P32));
P33(t) = vpa(subs(P33));

% NOTE: When using these expressions, when you plug in a value of t, the
% expression will still be a of class 'sym'.
% To evaluate these expressions as doubles use the following:
%       If t = 0:
%       P00_t1 = double(P00(0))
% Now P00_t1 will be a double

%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'Calculate the conditional probabilities as a function of rates {kij}';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);

%%
%     CALCULATE TCF's


% Need:
%       - Pji(t) calculated above
%       - Pi_eq equilibruim populations (from TDP_probcalc)
%           -> do not need TDP_probcalc - can calc Peq from condProb's
%       - Ai values of FRET states

t = sym('t');   % This allows for a dynamic workspace

% This is a vector for the FRET values - assign values
A1 = 0.42;              % These are just guesses for now
A2 = 0.51;
A3 = 0.67;

A = [ A1; A2; A3];

% Matrix of conditional probabilities Pi-->j with i is initial condition
cP = [ P11(t), P21(t), P31(t);
    P12(t), P22(t), P32(t);
    P13(t), P23(t), P33(t)];
% Row = final state & Column = initial condition

%%

% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P1EQ = P11(inf);
P2EQ = P22(inf);
P3EQ = P33(inf);

Peq = [P1EQ; P2EQ; P3EQ];

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


%--------------------------------------------------------------------------
% Calculate Two point TCF:
%-------------------------------------------------------------------------
C2sym(t) = 0*t;
for i = 1:numel(A)
    for j = 1:numel(A)
        %
        C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
        C2sym(t) = C2sym(t) + C2temp(t);
    end
end


if double(C2sym(0)) == double(msq)
    disp('Mean of the square <A^2> matches C2(t=0)!')
else
    disp('Problem: mean of the square <A^2> DOES NOT match C2(t=0)')
end

if double(C2sym(10^20)) == double(sqm)
    disp('Square of the mean <A>^2 matches C2(t=inf)!')
else
    disp('Problem: square of the mean <A>^2 DOES NOT match C2(t=inf)')
end

%--------------------------------------------------------------------------
% Evaluate C2 over a range of t's
%--------------------------------------------------------------------------

Npts = 150;
timeArray = [1:9,logspace(1,6.4771212,Npts)];
C2 = C2sym(timeArray);

%--------------------------------------------------------------------------
% Plot two point TCF
%--------------------------------------------------------------------------

%close all
figure(1)

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

saveName = ['C2_','example'];
saveas(TCF2pt,saveName, 'png')


%%
%--------------------------------------------------------------------------
% Four point TCF:
%--------------------------------------------------------------------------

% syms A0 A1 A2 A3
% sym t
t = sym('t');

% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P1EQ = double(P11(inf));
P2EQ = double(P22(inf));
P3EQ = double(P33(inf));

% global Peq
Peq = [P1EQ; P2EQ; P3EQ];

if sum(Peq) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('Problem: Equilibrium probabilities DO NOT sum to 1.')
end

% This is a vector for the FRET values - assign values
A1 = 0.42;% These are just guesses for now
A2 = 0.51;
A3 = 0.67;

A = [double(A1); double(A2); double(A3)]; % Double just makes sure these are numbers, and A is no longer symbolic.

Amean = sum(A.*Peq);
A = A - Amean;

msq = sum((A.^2).*Peq);     % square of mean <A^2>
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2


% Set tau 2
tau2 = 1;

% Create a vector of tau1 and tau3 values
Npts = 50;
tau1vec = logspace(0,3,Npts); % logspace(starting exponent, final exponent, number of bins)
tau3vec = logspace(0,3,Npts);



% Or define tau1 and tau3 symbolically and plot the functions
% tau1 = sym('tau1');
% tau2 = sym('tau2');

t1 = tau1vec;
t2 = tau2;
t3 = tau3vec';


tic
disp('... Calculating the Conditional Probabilities');
% Initialize size of 3D array for conditional probability for t1:
cP_t1 = zeros(numel(A),numel(A),length(t1));

% Define the vector in each position
cP_t1(1,1,:) = double(P11(t1));
cP_t1(2,1,:) = double(P12(t1));
cP_t1(3,1,:) = double(P13(t1));

cP_t1(1,2,:) = double(P21(t1));
cP_t1(2,2,:) = double(P22(t1));
cP_t1(3,2,:) = double(P23(t1));

cP_t1(1,3,:) = double(P31(t1));
cP_t1(2,3,:) = double(P32(t1));
cP_t1(3,3,:) = double(P33(t1));

% Initialize size of 3D array for conditional probability for t2:
cP_t2 = zeros(numel(A),numel(A),length(t2));

% Define the vector in each position ******* *****
cP_t2(1,1,:) = double(P11(t2));
cP_t2(2,1,:) = double(P12(t2));
cP_t2(3,1,:) = double(P13(t2));

cP_t2(1,2,:) = double(P21(t2));
cP_t2(2,2,:) = double(P22(t2));
cP_t2(3,2,:) = double(P23(t2));

cP_t2(1,3,:) = double(P31(t2));
cP_t2(2,3,:) = double(P32(t2));
cP_t2(3,3,:) = double(P33(t2));

% Initialize size of 3D array for conditional probability for t3:
cP_t3 = zeros(numel(A),numel(A),length(t3));

% Define the vector in each position
cP_t3(1,1,:) = double(P11(t3));
cP_t3(2,1,:) = double(P12(t3));
cP_t3(3,1,:) = double(P13(t3));

cP_t3(1,2,:) = double(P21(t3));
cP_t3(2,2,:) = double(P22(t3));
cP_t3(3,2,:) = double(P23(t3));

cP_t3(1,3,:) = double(P31(t3));
cP_t3(2,3,:) = double(P32(t3));
cP_t3(3,3,:) = double(P33(t3));


%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'calculate Conditional Probabilities.';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);


disp('... Calculating the 4 point TCF');
tic

% Create a matrix to hold the C4's calculated for each (tau1,tau3) pair for
% a set tau 2
C4 = zeros(length(tau1vec),length(tau3vec));
%-------------------------------------------------------------------------
% Iterate over all the Permutations of FRET States
%-------------------------------------------------------------------------
for i = 1:numel(A)
    for j = 1:numel(A)
        for k = 1:numel(A)
            for l = 1:numel(A)
                C4term_val =  A(l) *squeeze(cP_t3(l,k,:)) * A(k) * cP_t2(k,j) * A(j) * squeeze(cP_t1(j,i,:))'* A(i) * Peq(i);
                C4 = C4 + C4term_val;
            end
        end
    end
end

%Display the amount of time a process took. Begins at the last tic.
elapsedTime = toc;
task_str = 'calculate the four point TCF (C4).';
disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);


% Plot surface of C4
figure(2)
TCF4point = surf(tau1vec, tau3vec, C4);
title('Analytical Four-point TCF: C^{(4)}','FontSize',18)
xlabel('Time (\tau_1)','FontSize',14);
ylabel('Time (\tau_3)','FontSize',14);
zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';




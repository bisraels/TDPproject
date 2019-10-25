%--------------------------------------------------------------------------
% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: September 2019 (C4Maker_3state123_linear.m)
%
% PURPOSE:  Evaluate the linear 3-State (123) conditional probabilties with a set of
%           rates then the  2 point TCF and 4 point TCF for that system
%
% INPUT: (1) conditional probabilities from ODE solver: symCondProb_3state123_linear.mat
%
% OUTPUT: (1) FRET Histogram
%         (2) Two-point TCF C2
%         (3) Four-point TCF C4
%
% MODIFICATIONS:
%   (1) Adapted from C4Maker_3state123_cyclical.m
%
%----------------------------------------------------
function [C4,C4_diff,C2] = C4Maker_3state123_linear(t12,t21,t23,t32,A1,A2,A3,tau2,tau1range)
%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
verboseMode = 1; %Set to 1 to see alot of progress updates and print off.
clockMode = 0;
saveMode = 0;
plotMode = 1;

%MODEL: 1 <--> 2 <--> 3
programName = 'C4Maker_3state123_linear.m';
switch nargin
    case 0
        disp(['Using default values in ' programName]);
        
        t12_bounds = [1e-6,1000e-6];  %Paramater #1 is high--> med
        t21_bounds = [1e-6,1e-3];%Paramater #3 is med --> high
        t23_bounds = [1e-6,10e-3];%Paramater #4 is med --> low
        t32_bounds = [10e-6,10e-3];  %Paramater #5 is low --> Medium
        
        A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
        A2_bounds = [0.45,0.65];%Paramater #7 % Med FRET state
        A3_bounds = [0.30,0.45];%Paramater #8 %Low FRET state
        boundsArray = [t12_bounds;t21_bounds;t23_bounds;t32_bounds;A1_bounds;A2_bounds;A3_bounds];
        
        Nparams = length(boundsArray);
        population = rand(1,Nparams);
        for param_idx = 1:Nparams
            %To pick a random number in the interval of LB to UB:
            % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
            population(param_idx) = boundsArray(param_idx) + population(param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
        end
        t12 = population(1);
        t21 = population(2);
        t23 = population(3);
        t32 = population(4);
        A1 = population(5);
        A2 = population(6);
        A3 = population(7);
        
        tau2 = 0;
        
        Npts = 150;
        tau1range = [0:9,logspace(1,6.4771212,Npts)]/1e6;
    case 7
        
        tau2 = 0;
        
        Npts = 150;
        tau1range = [0:9,logspace(1,6.4771212,Npts)]/1e6;
end

%--------------------------------------------------------------------------
% Define the FRET Array
%--------------------------------------------------------------------------
A = [ A1; A2; A3];

%--------------------------------------------------------------------------
% Set the rates (4rates)
%--------------------------------------------------------------------------
k12 = 1/t12;
k21 = 1/t21;
k23 = 1/t23;
k32 = 1/t32;

% if verboseMode == 1
%     fprintf(['k12 = %f, k21 = %f, k23 = %f, k32 = %f'...
%         '\n A1 = %f, A2 = %f, A3 = %f\r\n'],...
%         k12,k21,k23,k32,A1,A2,A3);
% end
% %--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------

%Output created by ODE solver
if verboseMode == 1
    disp('Loading the conditional Probabilities as a function of rates');
end
if clockMode == 1, tic; end
load('symCondProb_3state123_linear.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')
% load('symCondProb_3state123_linear_soln1.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3');
%Display the amount of time a process      Took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'load the conditional probabilities';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

% disp('Calculating conditional probabilities using the rates defined')
t = sym('t');

%--------------------------------------------------------------------------
% CALCULATE THE EIGENVALUES
%--------------------------------------------------------------------------

if clockMode == 1
    tic;
end
% Evaluate the eigenvalues in terms of the rates defined above - produce as doubles
%subs(s) returns a copy of s, replacing symbolic variables in s, with their
%values obtained from the calling function and the MATLAB® Workspace,
% and then evaluates s. Variables with no assigned values remain as variables.
eval1 = double(vpa(subs(eval1)));
eval2 = double(vpa(subs(eval2)));
eval3 = double(vpa(subs(eval3)));

% if verboseMode == 1
%     disp(['The final eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/eval1) ' microseconds '...
%         ' tau2 = 1/eval2 = ' num2str(1e6*1/eval2) ' microseconds ' ...
%         'and tau3 = 1/eval3 = ' num2str(1e6*1/eval3) ' microseconds .']);
% end

if verboseMode == 1
    disp(['   ***The final eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/eval2) ' microseconds '...
        ' tau2 = 1/eval2 = ' num2str(1e6*1/eval3) ' microseconds.']);
end

%Display the amount of time a process      Took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'Calculate the eigenvalues';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end


%--------------------------------------------------------------------------
% Evaluate the conditional probabilities
%--------------------------------------------------------------------------
if clockMode == 1, tic; end
% Evaluate conditional probabilties by substituting in values from above
% and using vpa() to force the simplest form of the output.

%Pij(t) is prob from i--> j, assuming you start in state i: Pi(t=0)=100%=1

% P11(t) = vpa(subs(P11));
% P12(t) = vpa(subs(P12));
% P13(t) = vpa(subs(P13));
%
% P21(t) = vpa(subs(P21));
% P22(t) = vpa(subs(P22));
% P23(t) = vpa(subs(P23));
%
% P31(t) = vpa(subs(P31));
% P32(t) = vpa(subs(P32));
% P33(t) = vpa(subs(P33));

P11(t) = subs(P11);
P12(t) = subs(P12);
P13(t) = subs(P13);

P21(t) = subs(P21);
P22(t) = subs(P22);
P23(t) = subs(P23);

P31(t) = subs(P31);
P32(t) = subs(P32);
P33(t) = subs(P33);

%Display the amount of time a process      Took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'evaluate the conditional probabilities as a function of rates {kij}';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end
%
%--------------------------------------------------------------------------
% CALCULATE Equilibrium populations
%--------------------------------------------------------------------------

t = sym('t');   % This allows for a dynamic workspace

% Matrix of conditional probabilities Pi-->j with i is initial condition
cP = [ P11(t), P21(t), P31(t);
    P12(t), P22(t), P32(t);
    P13(t), P23(t), P33(t)];
% Row = final state & Column = initial condition
% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P1EQ = P11(inf);
P2EQ = P22(inf);
P3EQ = P33(inf);

Peq = [P1EQ; P2EQ; P3EQ];

%--------------------------------------------------------------------------
% (2) Calculate Two point TCF:
%-------------------------------------------------------------------------
%SUBTRACT OFF THE MEAN VALUE
Amean = sum(A.*Peq);
A = A - Amean;


if clockMode == 1, tic; end

C2sym(t) = 0*t;
for i = 1:numel(A)
    for j = 1:numel(A)
        %
        C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
        C2sym(t) = C2sym(t) + C2temp(t);
    end
end
%Display the amount of time a process      Took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'Calculate the 2-point TCF as a function of rates {kij}';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

msq = sum((A.^2).*Peq);     % square of mean <A^2>.
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2
if verboseMode == 1
    if double(C2sym(0)) == double(msq)
        disp('Mean of the square <A^2> matches C2(t=0)!')
    else
        fprintf('Problem: mean of the square <A^2> (%f) DOES NOT match C2(t=0) (%f)',double(msq),double(C2sym(0)))
    end
    
    if double(C2sym(10^20)) == double(sqm)
        disp('Square of the mean <A>^2 matches C2(t=inf)!')
    else
        
        fprintf('Problem: square of the mean <A>^2 (%f) DOES NOT match C2(t=inf) (%f)',double(sqm),double(C2sym(10^20)))
    end
end

%--------------------------------------------------------------------------
% Evaluate C2 over a range of t's
%--------------------------------------------------------------------------
C2 = C2sym(tau1range);
%--------------------------------------------------------------------------
%  Plot two point TCF
%--------------------------------------------------------------------------
%close all
if plotMode == 1
    figure(2)
    
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    %subplot(1,2,1)
    %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
    TCF2pt = plot(tau1range,C2,'LineWidth',2);
    
    title('Analytical Two point TCF','FontSize',18)
    xlabel('Time (\tau_1)','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);
    %xlim([10^-5 10^0])
    % ylim([C2(500) C2(0)])
    
    ax = gca;
    ax.XScale = 'log';
    set(gca,'yscale','log')
    
    if saveMode == 1
        saveName = ['C2_','example'];
        saveas(TCF2pt,saveName, 'png')
    end
end


disp('Will calculate C4');
tic
%--------------------------------------------------------------------------
% (3) Four point TCF:
%--------------------------------------------------------------------------

% Set tau 2

% Create a vector of tau1 and tau3 values
% Npts = 50;
% tau1vec = logspace(0,3,Npts); % logspace(starting exponent, final exponent, number of bins)
% tau3vec = logspace(0,3,Npts);
tau1vec = tau1range;
tau3vec = tau1range;

% Or define tau1 and tau3 symbolically and plot the functions
% tau1 = sym('tau1');
% tau2 = sym('tau2');

t1 = tau1vec;
t2 = tau2;
t3 = tau3vec';

if clockMode == 1, tic; end
if verboseMode == 1
    disp('... Calculating the Conditional Probabilities');
end
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

%Display the amount of time a process           Took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'calculate Conditional Probabilities.';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

if verboseMode
    disp(['... Calculating the 4 point TCF']);
end


% Create a matrix to hold the C4's calculated for each (tau1,tau3) pair for
% a set tau 2
if clockMode == 1
    tic;
end
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
C4 = double(C4);
%Display the amount of time a process           Took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'calculate the four point TCF (C4).';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

% disp(['Size C4 = ' num2str(size(C4))]);

C2Product = C2.*C2';
% disp(['Size C2Product = ' num2str(size(C2Product))]);
C4_diff = C4 - C2Product;
% disp(['Size C4_diff = ' num2str(size(C4_diff))]);

% Plot surface of C4
if plotMode == 1
    figure(3);
    
    set(gcf,'Color','w');
    set(gcf,'Name','C4');
    
    plot_TCF4point = surf(t1, t3, C4);
    title('Analytical Four-point TCF: C^{(4)}','FontSize',18)
    xlabel('Time (\tau_1)','FontSize',14);
    ylabel('Time (\tau_3)','FontSize',14);
    zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
    
    view(28,36);
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
end



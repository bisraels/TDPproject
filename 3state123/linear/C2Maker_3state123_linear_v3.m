%--------------------------------------------------------------------------
% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: September 2019 (C2Maker_3state123_linear.m)
%
% PURPOSE:  Evaluate the linear 3-State (123) conditional probabilties with a set of
%           rates then the  2 point TCF and 4 point TCF for that system
%
% INPUT: (1) conditional probabilities from ODE solver: symCondProb_3state123_linear.mat
%          (2) set of rates and fret states.
% OUTPUT: (1) C2 simulation
%--------------------------------------------------------------------------

% function C2_sim = C2Maker_3state123_linear_v3(t12,t21,t23,t32,A1,A2,A3,time)
%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
verboseMode = 1; %Set to 1 to see alot of progress updates and print off.
clockMode = 1;
saveMode = 0;
plotMode = 1;

useNumericalMethods_EqPop_mode = 0;


%MODEL: 1 <--> 2 <--> 3
programName = 'C2Maker_3state123_linear_v3.m';
% switch nargin
%     case 0
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

Npts = 150;
time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
%     case 7
%         Npts = 150;
%         time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
%
% end
%--------------------------------------------------------------------------
% Set the rates (4rates)
%--------------------------------------------------------------------------
k12 = 1/t12;
k21 = 1/t21;
k23 = 1/t23;
k32 = 1/t32;
K = [-k12, k21, 0;...
    k12, (-k21 - k23 ), k32;...
    0, k23, -k32;];
if verboseMode == 1
    fprintf(['k12 = %f, k21 = %f, k23 = %f, k32 = %f'...
        '\n A1 = %f, A2 = %f, A3 = %f\r\n'],...
        k12,k21,k23,k32,A1,A2,A3);
end

%--------------------------------------------------------------------------
% Define the FRET Array
%--------------------------------------------------------------------------
A = [ A1; A2; A3];

%--------------------------------------------------------------------------
% CALCULATE THE EIGENVECTORS AND EIGENVALUES
%--------------------------------------------------------------------------
if clockMode == 1
    tic;
end
[V,D] = eig(K);
%We needed to sort the eigenvalues because they are not in order
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

eval1 = Ds(3,3);        % First Eigenvalue should always be zero.
eval2 = Ds(2,2);
eval3 = Ds(1,1);


if verboseMode == 1
    disp(['The eigenvalues are: lam1 = ' num2str(eval1) ' sec^-1,'...
        'lam2 = ' num2str(eval2) ' sec^-1, '...
        'and lam3 = ' num2str(eval3) ' sec^-1.']);
end

if verboseMode == 1
    disp(['The eigentimescales are: tau1 = 1/lam1 = ' num2str(-1e6*1/eval1) ' microseconds '...
        ' tau2 = 1/eval2 = ' num2str(-1e6*1/eval2) ' microseconds ' ...
        'and tau3 = 1/eval3 = ' num2str(-1e6*1/eval3) ' microseconds.']);
end
if clockMode == 1
    task_str = 'calculate eigentimes';
    disp(['Took ' num2str(toc) ' seconds to ' task_str]);
end
%--------------------------------------------------------------------------
% CALCULATE Equilibrium populations (Numerical Methods)
%--------------------------------------------------------------------------

% if useNumericalMethods_EqPop_mode == 1
if clockMode == 1
    tic;
end
syms P1(t) P2(t) P3(t)

% Make a column vector of the equilibrium Populations
P(t) = [P1(t); P2(t); P3(t)];

% Define the DE we want to solve
eqn = diff(P(t),t)== K * P(t);
%OUTPUT: eqn =
%  diff(P1(t), t) == k21*P2(t) - k12*P1(t)
%  diff(P2(t), t) == k12*P1(t) - P2(t)*(k21 + k23) + k32*P3(t)
%  diff(P3(t), t) == k23*P2(t) - k32*P3(t)

% Solve the equations and give an output
%Pij(t) is prob from i--> j, assuming you start in state i: Pi(t=0)=100%=1
% Condition 1: P(t) = [P1(t) = 1, P2(t) = 0, P3(t) = 0]
Psol_1 = dsolve(eqn,[P1(0) == 1 , P2(0) == 0 , P3(0) == 0]);
%vpa(x) uses variable-precision floating-point arithmetic (VPA)
%The solution to the eigenvector problem will be included by dsove
P11(t) = vpa(Psol_1.P1);
% P12(t) = vpa(Psol_1.P2);
% P13(t) = vpa(Psol_1.P3);

% Condition 2: P(t) = [P1(t) = 0, P2(t) = 1, P3(t) = 0]
Psol_2 = dsolve(eqn,[P1(0) == 0 , P2(0) == 1 , P3(0) == 0]);
% P21(t) = vpa(Psol_2.P1);
P22(t) = vpa(Psol_2.P2);
% P23(t) = vpa(Psol_2.P3);

% Condition 3: P(t) = [P1(t) = 0, P2(t) = 0, P3(t) = 1]
Psol_1 = dsolve(eqn,[P1(0) == 0 , P2(0) == 0 , P3(0) == 1]);
% P31(t) = vpa(Psol_1.P1);
% P32(t) = vpa(Psol_1.P2);
P33(t) = vpa(Psol_1.P3);

t = sym('t');   % This allows for a dynamic workspace

% % Matrix of conditional probabilities Pi-->j with i is initial condition
% cP = [ P11(t), P21(t), P31(t);
%     P12(t), P22(t), P32(t);
%     P13(t), P23(t), P33(t)];
% % Row = final state & Column = initial condition
% % Define equilibrium populations from the conditional probabiltiies at
% % infinite time.

P1EQ = P11(100);
P2EQ = P22(100);
P3EQ = P33(100);

Peq = [P1EQ; P2EQ; P3EQ];
if clockMode == 1
    task_str = 'calculate equilibrium populations with Numerical methods';
    disp(['Took ' num2str(toc) ' seconds to ' task_str]);
end
%--------------------------------------------------------------------------
% CALCULATE Equilibrium populations (Analytical method)
%--------------------------------------------------------------------------

if clockMode == 1
    tic;
end
% Define equilibrium values of populations (altered, agrees with numerical)
P2EQ = 1/(k23/k32 + 1 + k21/k12);
P3EQ = k23/k32*P2EQ;
P1EQ = k21/k12*P2EQ;

Peq2 = [P1EQ; P2EQ; P3EQ];
if clockMode == 1
    task_str = 'calculate equilibrium populations with analytical methods';
    disp(['Took ' num2str(toc) ' seconds to ' task_str]);
end



%--------------------------------------------------------------------------
% (2) Calculate Two point TCF:
%-------------------------------------------------------------------------
Amean = sum(A.*Peq);
A = A - Amean;




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
    TCF2pt = plot(time,C2_sim,'LineWidth',2);
    
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


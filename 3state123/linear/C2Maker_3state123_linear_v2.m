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

% function C2_sim2 = C2Maker_3state123_linear_v2(t12,t21,t23,t32,A1,A2,A3,time)
%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
verboseMode = 0; %Set to 1 to see alot of progress updates and print off.
clockMode = 1;
saveMode = 0;
plotMode = 1;

use_C2_slow_numerical_mode = 1;
use_C2_fast_numerical_mode = 1;
use_C2_fast_analytical_mode = 0;
sim_4point_mode = 1;

%MODEL: 1 <--> 2 <--> 3
programName = 'C2Maker_3state123_linear.m';
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
% User Prefrences
%--------------------------------------------------------------------------
if clockMode == 1
    tic
end
%Output created by ODE solver
if verboseMode == 1
    disp('Loading the conditional Probabilities as a function of rates');
end
load('symCondProb_3state123_linear.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')

%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'load the conditional probabilities';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

if clockMode == 1
    tic
end
% disp('Calculating conditional probabilities using the rates defined')
t = sym('t');

%--------------------------------------------------------------------------
% CALCULATE THE EIGENVALUES
%--------------------------------------------------------------------------


if clockMode == 1
    tic
end
% Evaluate the eigenvalues in terms of the rates defined above - produce as doubles
%subs(s) returns a copy of s, replacing symbolic variables in s, with their
%values obtained from the calling function and the MATLAB� Workspace,
% and then evaluates s. Variables with no assigned values remain as variables.
eval1 = double(vpa(subs(eval1)));
eval2 = double(vpa(subs(eval2)));
eval3 = double(vpa(subs(eval3)));

if verboseMode == 1
    disp(['The final eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/eval1) ' microseconds '...
        ' tau2 = 1/eval2 = ' num2str(1e6*1/eval2) ' microseconds ' ...
        'and tau3 = 1/eval3 = ' num2str(1e6*1/eval3) ' microseconds .']);
end

%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'Calculate the eigenvalues';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end


%--------------------------------------------------------------------------
% Evaluate the conditional probabilities
%--------------------------------------------------------------------------

if clockMode == 1
    tic
end
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

%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'evaluate the conditional probabilities as a function of rates {kij}';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
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

if clockMode == 1
    tic
end
% Row = final state & Column = initial condition
% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P1EQ = P11(inf);
P2EQ = P22(inf);
P3EQ = P33(inf);

Peq = [P1EQ; P2EQ; P3EQ];
if clockMode == 1
    elapsedTime = toc;
    task_str = 'Calculate the Equilibrium Populations {Peq}';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end
%--------------------------------------------------------------------------
%  Calculate C2       (FAST ANALYTICAL UNTESTED)
%--------------------------------------------------------------------------
if use_C2_fast_analytical_mode == 1
    if clockMode == 1
        tic
    end
    K = [-k12, k21, 0;...
        k12, (-k21 - k23 ), k32;...
        0, k23, -k32;];
    [Vec, lambda] = eig(K);
    eval1 = lambda(1,1);        % First Eigenvalue should always be zero.
    eval2 = lambda(2,2);
    eval3 = lambda(3,3);
    
    % Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
    % basis [p2,p1,p0], where p# is the population in state #.
    a = (c1 + c2)/k21 - (k12 + k21)/k21;
    b = k12/k21 - (c1 + c2)/k21;
    c = 1;
    
    % Define eigenvector components corresponding to lam2: v2 = [x,y,z].
    x = (c1 - c2)/k21 - (k12 + k21)/k21;
    y = k12/k21 - (c1 - c2)/k21;
    z = 1;
    
    % Define equilibrium values of populations
    p2_eq = 1/(k21/k12 + 1 + k23/k32);
    p1_eq = k21/k12*p2_eq;
    p3_eq = k23/k32*p2_eq;
    
    % Define constants needed for differential equation solutions, for various
    % initial values of the population. m_2 is the constant "m" if the
    % population of state 2 begins at 1 (others would then be zero).
    n_3 = (-p2_eq - b/a*(1 - p3_eq))/(y - b/a*x);
    n_2 = (1 - p2_eq + b/a*p3_eq)/(y - b/a*x);
    n_1 = (-p2_eq + b/a*p3_eq)/(y - b/a*x);
    m_3 = (1 - p3_eq - n_3*x)/a;
    m_2 = (-p3_eq - n_2*x)/a;
    m_1 = (-p3_eq - n_1*x)/a;
    
    % Define populations as a function of time for different initial values.
    % pj_i population j as a function of time, with the population initially in
    % state i.
    p3_3 = p3_eq + m_3*a*exp(-lam1*time) + n_3*x*exp(-lam2*time);
    p3_2 = p3_eq + m_2*a*exp(-lam1*time) + n_2*x*exp(-lam2*time);
    p3_1 = p3_eq + m_1*a*exp(-lam1*time) + n_1*x*exp(-lam2*time);
    p2_3 = p2_eq + m_3*b*exp(-lam1*time) + n_3*y*exp(-lam2*time);
    p2_2 = p2_eq + m_2*b*exp(-lam1*time) + n_2*y*exp(-lam2*time);
    p2_1 = p2_eq + m_1*b*exp(-lam1*time) + n_1*y*exp(-lam2*time);
    p1_3 = p1_eq + m_3*c*exp(-lam1*time) + n_3*z*exp(-lam2*time);
    p1_2 = p1_eq + m_2*c*exp(-lam1*time) + n_2*z*exp(-lam2*time);
    p1_1 = p1_eq + m_1*c*exp(-lam1*time) + n_1*z*exp(-lam2*time);
    
    % Subtract mean values
    Amean = p1_eq*A1 + p2_eq*A2 + p3_eq*A3;
    A1 = A1 - Amean;
    A2 = A2 - Amean;
    A3 = A3 - Amean;
    
    % Calculate TCF
    % Calculate C2_sim
    C2_sim = p3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
        p2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
        p1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
    
    if clockMode == 1
        elapsedTime = toc;
        task_str = 'Calculate C2 usinng fast numerical methods';
        disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
    end
end

%--------------------------------------------------------------------------
% (2) Calculate Two point TCF: (slow numerical way)
%-------------------------------------------------------------------------
if use_C2_slow_numerical_mode == 1
    Amean = sum(A.*Peq);
    Ams = A - Amean;
    
    if clockMode == 1
        tic
    end
    
    C2sym(t) = 0*t;
    for i = 1:numel(Ams)
        for j = 1:numel(Ams)
            C2temp(t) = A(j) * cP(j,i) * Ams(i) * Peq(i);
            C2sym(t) = C2sym(t) + C2temp(t);
        end
    end
    
    C2_sim2 = C2sym(time);
    
    %Display the amount of time a process took. Begins at the last tic.
    if clockMode == 1
        elapsedTime = toc;
        task_str = 'Calculate the 2-point TCF as a function of rates {kij}';
        disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
    end
    
    
    %--------------------------------------------------------------------------
    %  Plot two point TCF
    %--------------------------------------------------------------------------
    %close all
    if plotMode == 1
        figure(2)
        clf;
        hold on;
        set(gcf,'Color','w');
        set(gcf,'Name','C2');
        %subplot(1,2,1)
        %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
        TCF2pt = plot(time,C2_sim2,'k-','LineWidth',2,'DisplayName','C2 (Numerical)');
        
        title('Two point TCF','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('C^{(2)}(\tau)','FontSize',14);
        %xlim([10^-5 10^0])
        % ylim([C2(500) C2(0)])
        
        ax = gca;
        ax.XScale = 'log';
        set(gca,'yscale','log')
        legend('show');
        if saveMode == 1
            saveName = ['C2_','example'];
            saveas(TCF2pt,saveName, 'png')
        end
    end
    
    %--------------------------------------------------------------------------
    % Sanity Checks
    %--------------------------------------------------------------------------
    msq = sum((Ams.^2).*Peq);     % square of mean <A^2>.
    sqm = (sum(Ams.*Peq))^2;      % mean square value <A>^2
    
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
end
%--------------------------------------------------------------------------
% (2) Calculate Two point TCF: (fast numerical way)
%-------------------------------------------------------------------------
if use_C2_fast_numerical_mode == 1
    Amean = sum(A.*Peq);
    Ams = A - Amean;
    
    %Use arrays of anonymous functions to calculate the conditional
    %probabilities
    cp = cell(3,3); % initialize a cell array
    %cp{i,j} = @t cond. prob from i to j
    cp{1,1} = @(t,k12,k21,k23,k32) (k21*k32)/(k12*k23 + k12*k32 + k21*k32) - (0.5*k12*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*(k23 + k32))/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 + 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) + (0.5*k12*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*(k23 + k32))/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 - 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{1,2} = @(t,k12,k21,k23,k32) (k12*k32)/(k12*k23 + k12*k32 + k21*k32) + (0.5*k12*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*k32)/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 + 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) - (0.5*k12*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*k32)/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 - 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{1,3} = @(t,k12,k21,k23,k32) (k12*k23)/(k12*k23 + k12*k32 + k21*k32) + (0.5*k12*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) - (0.5*k12*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{2,1} = @(t,k12,k21,k23,k32) (k21*k32)/(k12*k23 + k12*k32 + k21*k32) + (0.5*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*(k23 + k32))/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 + 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12*k23 - 1.0*k12*k21 + k12*k32 + 2.0*k21*k32 + k12*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) - 1.0*k12^2))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) + (0.5*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*(k23 + k32))/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 - 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12*k21 - 1.0*k12*k23 - 1.0*k12*k32 - 2.0*k21*k32 + k12*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k12^2))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{2,2} = @(t,k12,k21,k23,k32) (k12*k32)/(k12*k23 + k12*k32 + k21*k32) - (0.5*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*k32)/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 - 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12*k21 - 1.0*k12*k23 - 1.0*k12*k32 - 2.0*k21*k32 + k12*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k12^2))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) - (0.5*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*k32)/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 + 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12*k23 - 1.0*k12*k21 + k12*k32 + 2.0*k21*k32 + k12*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) - 1.0*k12^2))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{2,3} = @(t,k12,k21,k23,k32) (k12*k23)/(k12*k23 + k12*k32 + k21*k32) - (0.5*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*(k12*k23 - 1.0*k12*k21 + k12*k32 + 2.0*k21*k32 + k12*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) - 1.0*k12^2))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) - (0.5*k23*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*(k12*k21 - 1.0*k12*k23 - 1.0*k12*k32 - 2.0*k21*k32 + k12*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k12^2))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{3,1} = @(t,k12,k21,k23,k32) (k21*k32)/(k12*k23 + k12*k32 + k21*k32) - (0.5*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*(k23 + k32))/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 - 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12^2*k32 - 1.0*k12*k32^2 - 1.0*k21*k32^2 + k21^2*k32 + 2.0*k12*k21*k32 - 1.0*k12*k23*k32 + k21*k23*k32 + k12*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k21*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) - (0.5*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*(k23 + k32))/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 + 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12*k32^2 - 1.0*k12^2*k32 + k21*k32^2 - 1.0*k21^2*k32 - 2.0*k12*k21*k32 + k12*k23*k32 - 1.0*k21*k23*k32 + k12*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k21*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{3,2} = @(t,k12,k21,k23,k32) (k12*k32)/(k12*k23 + k12*k32 + k21*k32) + (0.5*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*k32)/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 + 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12*k32^2 - 1.0*k12^2*k32 + k21*k32^2 - 1.0*k21^2*k32 - 2.0*k12*k21*k32 + k12*k23*k32 - 1.0*k21*k23*k32 + k12*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k21*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) + (0.5*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*((1.0*k32)/k23 - (0.5*k12 + 0.5*k21 + 0.5*k23 + 0.5*k32 - 0.5*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2))/k23)*(k12^2*k32 - 1.0*k12*k32^2 - 1.0*k21*k32^2 + k21^2*k32 + 2.0*k12*k21*k32 - 1.0*k12*k23*k32 + k21*k23*k32 + k12*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k21*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    cp{3,3} = @(t,k12,k21,k23,k32) (k12*k23)/(k12*k23 + k12*k32 + k21*k32) + (0.5*exp(-0.5*t*(k12 + k21 + k23 + k32 - 1.0*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*(k12^2*k32 - 1.0*k12*k32^2 - 1.0*k21*k32^2 + k21^2*k32 + 2.0*k12*k21*k32 - 1.0*k12*k23*k32 + k21*k23*k32 + k12*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k21*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)) + (0.5*exp(-0.5*t*(k12 + k21 + k23 + k32 + (k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))*(k12*k32^2 - 1.0*k12^2*k32 + k21*k32^2 - 1.0*k21^2*k32 - 2.0*k12*k21*k32 + k12*k23*k32 - 1.0*k21*k23*k32 + k12*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2) + k21*k32*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2)))/((k12*k23 + k12*k32 + k21*k32)*(k12^2 + 2.0*k12*k21 - 2.0*k12*k23 - 2.0*k12*k32 + k21^2 + 2.0*k21*k23 - 2.0*k21*k32 + k23^2 + 2.0*k23*k32 + k32^2)^(1/2));
    
    if clockMode == 1
        tic
    end
    t = time;
    C2_sim2 = zeros(1,length(time));
    for i = 1:numel(Ams)
        for j = 1:numel(Ams)
%             C2_sim_temp = Ams(j) .* cp{i,j}(time,k12,k21,k23,k32) * Ams(i) * Peq(i);
            C2_sim_temp = Ams(j) .* double(subs(cp{j,i})) * Ams(i) * Peq(i);
            C2_sim2 = C2_sim2 + C2_sim_temp;
        end
    end
    %Display the amount of time a process took. Begins at the last tic.
    if clockMode == 1
        elapsedTime = toc;
        task_str = 'Calculate the 2-point TCF using fast numerical methods';
        disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
    end
    
    
    
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
        TCF2pt = plot(time,C2_sim2,'r--','LineWidth',2,'DisplayName','C2 (fast Numerical)');
        
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
end

class(cp{1,1})

%--------------------------------------------------------------------------
%  Plot four point TCF (Fast Numerical Working)
%--------------------------------------------------------------------------
if sim_4point_mode == 1
    C4_sim = zeros(numel(C2_sim2),numel(C2_sim2));
    C4term_val = zeros(size(C4_sim));
    tau2 = 0;
    tic
    tocArray = zeros(1,numel(Ams)^4);
    for i = 1:numel(Ams)
        for j = 1:numel(Ams)
            for k = 1:numel(Ams)
                for l = 1:numel(Ams)
                    C4term_val =  Ams(l) .* cp{k,l}(time,k12,k21,k23,k32) * Ams(k) .* cp{j,k}(tau2,k12,k21,k23,k32) * Ams(j) .* cp{i,j}(time,k12,k21,k23,k32)'* Ams(i) * Peq(i);
                    C4_sim = C4_sim + C4term_val;
                    num = 3^3*(i-1)+3^2*(j-1)+3^1*(k-1)+3^0*(l-1)+1;
                    disp(num2str(num));
                    tocArray(num) =  toc;
                end
            end
        end
    end
    C4_sim = double(C4_sim);
    %Display the amount of time a process           Took. Begins at the last tic.
    if clockMode == 1
        elapsedTime = toc;
        task_str = 'calculate the four point TCF (C4).';
        disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
    end
    if plotMode == 1
        figure(4);
        plot([1:numel(tocArray)],tocArray);
        xlabel('Number of calcs');
        ylabel('time');
    end
    
    % disp(['Size C4 = ' num2str(size(C4))]);
    
    C2_sim_Product = C2_sim2.*C2_sim2';
    % disp(['Size C2Product = ' num2str(size(C2Product))]);
    C4_diff_sim = C4_sim - C2_sim_Product;
    % disp(['Size C4_diff = ' num2str(size(C4_diff))]);
    
    % Plot surface of C4
    if plotMode == 1
        figure(3);
        hold on;
        set(gcf,'Color','w');
        set(gcf,'Name','C4');
        
        plot_TCF4point = mesh(time, time, C4_sim);
        title('Analytical Four-point TCF: C^{(4)}','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('Time (\tau_3)','FontSize',14);
        zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
        
        view(28,36);
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
    end
end
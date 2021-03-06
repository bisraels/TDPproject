%--------------------------------------------------------------------------
% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: September 2019 (histMaker_3state123_cyclical.m)
%
% PURPOSE:  Evaluate the Linear 3-State (123) conditional probabilties with a set of
%           rates then the  2 point TCF and 4 point TCF for that system
%
% INPUT: (1) conditional probabilities from ODE solver: symCondProb_3state123_cyclical.mat
%
% OUTPUT: (1) FRET Histogram
%         (2) Two-point TCF C2
%         (3) Four-point TCF C4
%
% MODIFICATIONS:
%   (1) Adapted from function_3state123_cyclical.m
%
%--------------------------------------------------------------------------

function [Peq] = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3)
programName = 'histMaker_3state123_cyclical';
%--------------------------------------------------------------------------
% SET PARAMATERS
%--------------------------------------------------------------------------
switch nargin
    case 0
        disp(['Using default values in ' programName]);
        
        t12_bounds = [1e-6,1000e-6];  %Paramater #1 is high--> med
        t13_bounds = [100e-6,10e-3];    %Paramater #2 is high --> low
        t21_bounds = [1e-6,1e-3];%Paramater #3 is med --> high
        t23_bounds = [1e-6,10e-3];%Paramater #4 is med --> low
        t31_bounds = [10e-6,10e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
        % *t32 wll be determined by the other rates
        
        A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
        A2_bounds = [0.45,0.65];%Paramater #7 % Med FRET state
        A3_bounds = [0.30,0.45];%Paramater #8 %Low FRET state
        boundsArray = [t12_bounds;t13_bounds;t21_bounds;t23_bounds;t31_bounds;A1_bounds;A2_bounds;A3_bounds];
         
        Nparams = length(boundsArray);
        population = rand(1,Nparams);
        for param_idx = 1:Nparams
            %To pick a random number in the interval of LB to UB:
            % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
            population(param_idx) = boundsArray(param_idx) + population(param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
        end
        t12 = population(1);
        t13 = population(2);
        t21 = population(3);
        t23 = population(4);
        t31 = population(5);
        A1 = population(6);
        A2 = population(7);
        A3 = population(8);
end

%--------------------------------------------------------------------------
% Define the FRET Array
%--------------------------------------------------------------------------

A = [ A1; A2; A3];

%--------------------------------------------------------------------------
% Set the rates
%--------------------------------------------------------------------------
k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;
% % Detailed balance condition: %k31 will be the rate fixed by the others
k32 = k12*k23*k31/(k13*k21);
% t32 = 1/k32;
%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
verboseMode = 0; %Set to 1 to see alot of progress updates and print off.
clockMode = 0;
saveMode = 0;
plotMode = 0;

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
% load('symCondProb_3state123_cyclical.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')
load('symCondProb_3state123_cyclical.mat','P11','P22','P33')

%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'load the conditional probabilities';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end



% 
% %--------------------------------------------------------------------------
% % CALCULATE THE EIGENVALUES
% %--------------------------------------------------------------------------
% 
% % Evaluate the eigenvalues in terms of the rates defined above - produce as doubles
% %subs(s) returns a copy of s, replacing symbolic variables in s, with their
% %values obtained from the calling function and the MATLAB� Workspace,
% % and then evaluates s. Variables with no assigned values remain as variables.
% eval1 = double(vpa(subs(eval1)));
% eval2 = double(vpa(subs(eval2)));
% eval3 = double(vpa(subs(eval3)));
% 
% if verboseMode == 1
%     disp(['The final eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/eval1) ' microseconds '...
%         ' tau2 = 1/eval2 = ' num2str(1e6*1/eval2) ' microseconds ' ...
%         'and tau3 = 1/eval3 = ' num2str(1e6*1/eval3) ' microseconds .']);
% end
% %Display the amount of time a process took. Begins at the last tic.
% if clockMode == 1
%     elapsedTime = toc;
%     task_str = 'calculate the eigenvalues';
%     disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
% end
% 
% % Evaluate conditional probabilties by substituting in values from above
% % and using vpa() to force the simplest form of the output.
% 
% %Pij(t) is prob from i--> j, assuming you start in state i: Pi(t=0)=100%=1
if clockMode == 1
tic
end
% disp('Calculating conditional probabilities using the rates defined')
t = sym('t');
P11(t) = vpa(subs(P11));
% P12(t) = vpa(subs(P12));
% P13(t) = vpa(subs(P13));
% 
% P21(t) = vpa(subs(P21));
P22(t) = vpa(subs(P22));
% P23(t) = vpa(subs(P23));
% 
% P31(t) = vpa(subs(P31));
% P32(t) = vpa(subs(P32));
P33(t) = vpa(subs(P33));
% 
% % NOTE: When using these expressions, when you plug in a value of t, the
% % expression will still be a of class 'sym'.
% % To evaluate these expressions as doubles use the following:
% %       If t = 0:
% %       P00_t1 = double(P00(0))
% % Now P00_t1 will be a double
% 
% %Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'Calculate the conditional probabilities as a function of rates {kij}';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end
% %
% %     CALCULATE TCF's
% 
% 
% % Need:
% %       - Pji(t) calculated above
% %       - Pi_eq equilibruim populations (from TDP_probcalc)
% %           -> do not need TDP_probcalc - can calc Peq from condProb's
% %       - Ai values of FRET states
% 
% t = sym('t');   % This allows for a dynamic workspace
% 
% % Matrix of conditional probabilities Pi-->j with i is initial condition
% cP = [ P11(t), P21(t), P31(t);
%     P12(t), P22(t), P32(t);
%     P13(t), P23(t), P33(t)];
% % Row = final state & Column = initial condition

%
if clockMode == 1
tic
end
% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P1EQ = P11(inf);
P2EQ = P22(inf);
P3EQ = P33(inf);
%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'calculate the equilibrium populations by evaluating Pii(inf)';
    disp(['     Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

Peq = [P1EQ; P2EQ; P3EQ];

if verboseMode == 1
    if sum(Peq) == 1
        disp('Equilibrium probabilities sum to 1!')
    else
        error('Problem: Equilibrium probabilities DO NOT sum to 1.')
    end
end

if plotMode == 1
    figure(1);
    
    set(gcf,'Color','w');
    set(gcf,'Name','FRET Histogram');
%     set(gcf,'Position',[1681 437 631 511]);
%     
    sigma_A1 = 0.15;
    sigma_A2 = 0.1;
    sigma_A3 = 0.1;
    
    % Final fit histogram
    FRET_bins = linspace(0,1,100);
    hist_sim = P1EQ*exp(-((FRET_bins-A1)/sigma_A1).^2) + P2EQ*exp(-((FRET_bins-A2)/sigma_A2).^2) + P3EQ*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim= sum(hist_sim);
    
    hist_sim = hist_sim./denom_hist_sim;
    
    % One line for the final fit
    plot(FRET_bins,hist_sim,'r-','LineWidth',2,'DisplayName','Final Fit');
    xlabel('FRET Efficiency','FontSize',14);
    ylabel('Frequency','FontSize',14);
    title('Simulated Histograms','FontSize',14);
    
    
    hold on;
    % Histogram for each population individually
    plot(FRET_bins,P1EQ*exp(-((FRET_bins-A1)/sigma_A1).^2)./denom_hist_sim,'c--','LineWidth',1,'DisplayName','P1_{eq}');
    plot(FRET_bins,P2EQ*exp(-((FRET_bins-A2)/sigma_A2).^2)./denom_hist_sim,'m--','LineWidth',1,'DisplayName','P2_{eq}');
    plot(FRET_bins,P3EQ*exp(-((FRET_bins-A3)/sigma_A3).^2)./denom_hist_sim,'g--','LineWidth',1,'DisplayName','P3_{eq}');
    lgd = legend('show');
    lgd.Location = 'northwest';
    lgd.FontSize = 14;
    hold off;
    
    set(gca,'FontSize',14);
end


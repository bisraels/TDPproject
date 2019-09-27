%__________________________________________________________________________
% AUTHOR: Brett Israels (September 2019 Creation)
%          - Inspired by Carey Phelp's GenAlg_v16_0p5uMgp32.m =
%
% NAME: GenAlg_3state123_cyclical_exp.m
%
% FUNCTION: % Fits (1) the FRET Histogram, (2) the 2pt TCF and (3) the 4-point TCF
%
% PROCEDURE:
% (1) Load the experimental histogram
% (2) Load the experimental 2pt time correlation functions to fit
% (3) Load the experimental 4pt TCFs for various tau2 values
% (4) Make an initial population of guesses for the model paramaters
% (5) Construct Histograms, C2, and {C4} for each guess
% (6) Compare each simulated plot to the experimental one to calculate rms
% (7) Choose the best elements of each and mix together.
% (8) Repeat steps 5-7 until the simulation converges below the threshold
%
% INPUT:
% (1) FRET Histogram : '3p15mer_0p0uMgp32_GaussianFit_N3.mat'
% (2) Two-Point TCF : '3p15mer_0p0uMgp32_tcfavg_2exp_fitresult.mat'
% (3) Four-Point TCF : '000010us_tau2-000000us_Neq035_fourptTCFavg.mat'
%
% OUTPUT:
% (1) BestFitResults.mat         %All the fitting paramaters
% (2) BestFitRestults_hist.fig   %
% (3) fitInputData.mat
% (4) genAlgParamaters.mat
% (5) ModelResultsFigure.fig
% (6) plottingParamaters.mat
%
% CALLS:
% (1) ODEsolver_3state123_cyclical.m %Calculates the conditional probabilities
% (2) function_3state123_cyclical.m  % Calculates the Histogram/C2/C4
%
% (1) C2_sim = TCF_cyclic3state                                 % Calculates C2
% (2) [p1_eq,p2_eq,p3_eq] = cyclic3state_hist                   % Calculates the histogram values
% (3)[C4_tau2eq0_sim,~] = FourPtTCF_cyclic3state_norm           % Calculates the 4 point TCF
%
% MODIFICATION LOG:
% BI 20190924 Added Clair and I's Updated scripts to calculate Eq Pop.
% BI 20190926 Added the yoffset to the simulations of C2 (for better chisquared comparison)
% BI 20190926
%__________________________________________________________________________

%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
tic
% function guess = GenAlg_v17_0p0uMgp32
programName = 'GenAlg_3state123_cyclical_exp';
disp(['Now Running ' programName '.m']);
close all


%//////////////////////////////////////////////////////////////////////////
% PART 1: Set the program up with various options
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
global normalizeMode verboseMode guessUpdateMode diagnoseMode

normalizeMode = 0;
verboseMode = 0;
guessUpdateMode = 0;
diagnoseMode = 1;

clockMode = 0;
saveMode = 0;

plotMode = 1;
plot_Gen1ChiSquaredMode = 1;%(4)  Adds a point to the chisquare vs member plot(plot4) SLOW
plot_Gen1_memberGuessesMode = 0; % (1-2-3) Plots Histograms, C2, C4 of gen1 (SLOW)
plot_GenerationalProgressMode = 1;%(4) Gray circles on figure(4)
plot_ChiSquaredHistogramMode = 1;%(7) Histogram of chisquare values
pauseBetweenGenerationMode = 0;
useFigurePosnMode = 1;
showBestFitMode = 1;

global fitHistMode fitC2Mode fitC4Mode useOldCodeMode
fitHistMode = 1; fitHistData_mode = 0;%If 0 you will fit to the histogram fit
fitC2Mode = 1;
fitC4Mode = 0;
if sum([fitHistMode,fitC2Mode,fitC4Mode]) ==0
    error('No surfaces to optimize to! Pick atleast one.');
end

%**************************************************************************
useOldCodeMode = 	0;
if useOldCodeMode == 1
    disp('     ***USING OLD CODES');
else
    disp('     ***Using NEW codes');
    % if verboseMode == 1
    disp('Loading the conditional Probabilities as a function of rates');
    % end
    %Output created by ODE solver. Each variable is a function of rates {kij}
    load('symCondProb_3state123_cyclical.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')
end
%**************************************************************************
%--------------------------------------------------------------------------
% Declare global variables
%--------------------------------------------------------------------------
global genNum targetHistogram weightingFactor_FREThist FRET_bins
global C2_exp_x C2_exp_y weightingFactor_C2  weightC2func
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func

%Choose Weighting Amount
weightingFactor_FREThist = 10;%10*3*300; % Weighting for FRET hist comparison
weightingFactor_C2 = 1;%10*3*10000;
weightingFactor_C4_t0 = 1;%1*5000;

%--------------------------------------------------------------------------
% Genetic Algorithm Patamaters
%--------------------------------------------------------------------------
% Fitting real expt0.
NmembersInitPop = 100;
maxGenerations = 200;
maxRepeats = 10;

percentReproduce = .2;
Nreproduce = round(percentReproduce*NmembersInitPop);    % Number of individuals who reproduce per generation.
fprintf('     Each Generation will have %d/%d of its members reproduce\r',Nreproduce,NmembersInitPop);

Nmutations = 1*NmembersInitPop;                   % Number of mutations per generation.
fprintf('     Each Generation will see a total of %d mutations. (%d per individual on average).\r',Nmutations,Nmutations/NmembersInitPop);

percentToKeep = 2;
membersToKeep = NmembersInitPop*percentToKeep/100;
fprintf('     Each Generation will have the top %d percent unmolested (%d/%d individuals).\r',percentToKeep,membersToKeep,NmembersInitPop);

threshold = 0.001;                          % Minimal allowable percentage difference before program quits.


%//////////////////////////////////////////////////////////////////////////
% PART 1: Model Specific Paramaters
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% CYCLICAL THREE STATE MODEL : t32 is fixed by the other rates.
% 1 <--> 2 <--> 3 <--> 1
global sigma_A1 sigma_A2 sigma_A3 %These will be used in the optimization
sigma_A1 = 0.15;
sigma_A2 = 0.1;
sigma_A3 = 0.1;

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

Nparam = length(boundsArray);
maxMutationCountsArray = zeros(1,Nparam);
minMutationCountsArray = zeros(1,Nparam);

% Maximum mutation value
Max_mut_factor = .1;
t12_mutate = Max_mut_factor*t12_bounds(2)/2;%100e-6;
t13_mutate = Max_mut_factor*t13_bounds(2)/2;%100e-6;
t21_mutate = Max_mut_factor*t21_bounds(2)/2;%500e-6;
t23_mutate = Max_mut_factor*t23_bounds(2)/2;%10e-6;
t31_mutate = Max_mut_factor*t31_bounds(2)/2;%100e-6;
A1_mutate = (A1_bounds(2) - A1_bounds(1))*0.1;%0.05;
A2_mutate = (A2_bounds(2) - A2_bounds(1))*0.1;%0.05;0.05;
A3_mutate = (A3_bounds(2) - A3_bounds(1))*0.1;%0.05;0.05;
maxMutationArray = [t12_mutate;t13_mutate;t21_mutate;t23_mutate;t31_mutate;A1_mutate;A1_mutate;A2_mutate];

%//////////////////////////////////////////////////////////////////////////
% PART 2: Load the target data
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%Start in the single molecule folder (smData_Processed): comp specific
wd = pwd;
if contains(wd,'C:\Users\baimi\')%Work computer
    cd('C:\Users\baimi\Dropbox\MarcusLab\Data\smData_Processed\S1S2\gp32_0p0uM\ChosenMolecules\genAlgFiles');
    computer_baimi_mode = 1;
    computer_terminal_str = 'computer_baimi_mode';
    disp(['Turning on on ' computer_terminal_str]);
elseif contains(wd,'/Users/bisraels')%macbook
    cd('/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/S1S2/gp32_0p0uM/ChosenMolecules/genAlgFiles');
    computer_bisraels_mode = 1;
    computer_terminal_str = 'computer_bisraels_mode';
    disp(['Turning on on ' computer_terminal_str]);
else
    general_computer_mode = 1;
    computer_terminal_str = 'computer_general_mode';
    disp(['Turning on on ' computer_terminal_str]);
end

wd = pwd;
disp(['You are now in the realm of ' wd]);

%--------------------------------------------------------------------------
% Filenames of the data (PART 2: Load the target data)
%--------------------------------------------------------------------------
histogram_FileName = '3p15mer_0p0uMgp32_GaussianFit_N3.mat';
C2_FileName = '3p15mer_0p0uMgp32_tcfavg_2exp_fitresult.mat';
C4_FileName = '000010us_tau2-000000us_Neq035_fourptTCFavg.mat';

%--------------------------------------------------------------------------
% (1) Optimization Target #1: 1-D FRET HISTOGRAM (PART 2: Load the target data)
%--------------------------------------------------------------------------
if fitHistMode == 1
    load(histogram_FileName,'xData','yData','xFit','yFit');
    if fitHistData_mode == 1
        FRET_bins = xData;
        targetHistogram = yData;
    elseif fitHistData_mode == 0
        FRET_bins = xFit;
        targetHistogram = yFit;
    end
    
    targetHistogram = targetHistogram./sum(targetHistogram);
    
    if plotMode
        figure(1)
        set(gcf,'Name','FRET Histogram');
        if useFigurePosnMode == 1
            switch computer_terminal_str
                case 'computer_baimi_mode'
                    set(gcf,'Position',[1963 573 560 420]);%Work Desktop
                case 'computer_bisraels_mode'
                    set(gcf,'Position',[560 535 560 420]);%Macbook pro
            end
            
        end
        set(gcf,'Color','w');
        
        plot(FRET_bins,targetHistogram,'b-','DisplayName','Data');
        xlabel('FRET Efficiency','FontSize',14);
        ylabel('Frequency','FontSize',14);
        title('Experimental vs Simulated Histograms','FontSize',14);
        lgd = legend('show');
        lgd.Location = 'northwest';
        drawnow();
    end
end

%--------------------------------------------------------------------------
% (2) Optimization Target #2: 2-point TCF (C2)           (20 uSec and on)  (PART 2: Load the target data)
%--------------------------------------------------------------------------
if fitC2Mode == 1
    load(C2_FileName,'time','yData','y','yoff');
    fitC2Data_mode = 1;%If 0 you will fit to the histogram fit
    if fitC2Data_mode == 1
        C2_exp_x = time;
        C2_exp_y = yData;
    elseif fitC2Data_mode == 0
        C2_exp_x = time;
        C2_exp_y = y;
    end
    
    
    if normalizeMode == 1
        C2_exp_y = C2_exp_y./C2_exp_y(1);
    end
    
    % 2-pt. TCF weighting function
    weightC2func = 1./(sqrt(C2_exp_x));
    
    
    if plotMode == 1
        figure(2)
        clf;
        set(gcf,'Name','Two-Point TCF: C^{(2)}');
        if useFigurePosnMode == 1
            set(gcf,'Position',[2535 571 560 420]);
        end
        set(gcf,'Color','w');
        
        plot(C2_exp_x,C2_exp_y,'b.','MarkerSize',10,'DisplayName','C^{(2)}(\tau) Data');
        xlabel('Time (sec)','FontSize',14);
        ylabel('C2','FontSize',14);
        title('Experimental vs Simulated C2','FontSize',14);
        set(gca,'xscale','log');
        %         legend('show');
        drawnow();
    end
end
%--------------------------------------------------------------------------
% (3) Optimization Target #3: 4-point TCF (C4)           (10 uSec and on) (PART 2: Load the target data)
%--------------------------------------------------------------------------
if fitC4Mode == 1
    load(C4_FileName,'FourPtTCF_avg','tau1arrayUsec','tau3arrayUsec','tau2ValUsec');
    C4_tau2eq0_exp = FourPtTCF_avg;
    C4_tau1range = tau1arrayUsec*1e-6;
    C4_tau3range = tau3arrayUsec*1e-6';
    
    if normalizeMode == 1
        C4_tau2eq0_exp = C4_tau2eq0_exp./C4_tau2eq0_exp(1,1);
    end
    
    % 4-pt. TCF weighting function
    wC4func = 1./sqrt(C4_tau1range).*(1./(sqrt(C4_tau3range)));
    
    if plotMode == 1
        figure(3)
        clf;
        set(gcf,'Name','Four-Point TCF: C^{(4)}');
        if useFigurePosnMode == 1
            set(gcf,'Position', [3117 569 560 420]);
        end
        set(gcf,'Color','w');
        
        C4dataPlot = mesh(C4_tau1range,C4_tau3range,C4_tau2eq0_exp,'DisplayName','C^{(4)} Data');
        xlabel('\tau_1 (sec)','FontSize',14);
        ylabel('\tau_3 (sec)','FontSize',14);
        ylabel('C^{(4)}','FontSize',14);
        title('Experimental vs Simulated C4','FontSize',14);
        set(gca,'xscale','log');
        drawnow();
        
        
        [sample_description, save_prefix] = sample_descriptionGetter();
        title_str = ['C^{(4)}(\tau_1, \tau_2 = ' num2str(tau2ValUsec) '\musec, \tau_3)' ...
            10 sample_description];
        title(title_str,'fontsize',20);
        xlabel('\tau_1 (sec)','fontsize',18);
        ylabel('\tau_3 (sec)','fontsize',18);
        zlabel('C^{(4)}(\tau_1, \tau_2, \tau_3)','fontsize',18);
        
        %----------| Clean up the plot (non-specific --> specific) |-----------
        grid on;
        set(gca,'FontSize',12);
        axis tight;
        axis square;
        colormap jet;
        colorbar;
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        %         set(gca,'zscale','log');
        %         zlim([1e-4,inf]);
        view(55.3,33.2);
        
        drawnow();
    end
end

%//////////////////////////////////////////////////////////////////////////
% rmscalc will be the handle of the function which computes the ChiSquare
%//////////////////////////////////////////////////////////////////////////
% rmscalc = @multigoaltcf_analytical_3state;

%//////////////////////////////////////////////////////////////////////////
% PART 3: Make the initial generation of guesses
%//////////////////////////////////////////////////////////////////////////
% gen1
%Give each member of the population some value between the LB and UB
population = rand(NmembersInitPop,Nparam);
for param_idx = 1:Nparam
    %To pick a random number in the interval of LB to UB:
    % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
    population(:,param_idx) = boundsArray(param_idx,1) + population(:,param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
end


newpopulation = zeros(size(population));
pop_chisquared = zeros(NmembersInitPop,1);
guess = zeros(maxGenerations,9);
genNum = 1;
for pop_idx = 1:NmembersInitPop
    t12 = population(pop_idx,1);
    t13 = population(pop_idx,2);
    t21 = population(pop_idx,3);
    t23 = population(pop_idx,4);
    t31 = population(pop_idx,5);
    A1 = population(pop_idx,6);
    A2 = population(pop_idx,7);
    A3 = population(pop_idx,8);
    
    
    %--------------------------------------------------------------------------
    % gen1: Plot the array of initial guesses
    %--------------------------------------------------------------------------
    if plot_Gen1_memberGuessesMode == 1
        %--------------------------------------------------------------------------
        % (1) Optimization Target #1: 1-D FRET HISTOGRAM  (Plot the initial
        % histograms which would result from choices (PART 3: Make the % initial generation of guesses)
        %--------------------------------------------------------------------------
        if fitHistMode == 1
            figure(1)
            hold on;
            if useOldCodeMode == 1
                k12 = 1/t12;
                k13 = 1/t13;
                k21 = 1/t21;
                k23 = 1/t23;
                k31 = 1/t31;
                k32 = k12*k23*k31/(k13*k21);
                % %         function [p0_eq,p1_eq,p2_eq,lam1,lam2] = FourPtTCF_cyclic3state_xo(tau1range,tau2,A0,A1,A2,k01,k10,k12,k21,k20)
                %                 [p1_eq,p2_eq,p3_eq,~,~] = FourPtTCF_cyclic3state_xo(C4_tau1range,0,A1,A2,A3,k12,k21,k23,k32,k31);
                [p1_eq,p2_eq,p3_eq] = cyclic3state_hist(A1,A2,A3,k12,k21,k23,k32,k31);
            elseif useOldCodeMode == 0
                [Peq] = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
                p1_eq = Peq(1);
                p2_eq = Peq(2);
                p3_eq = Peq(3);
            end
            
            hist_sim = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
            denom_hist_sim = sum(hist_sim);
            hist_sim = hist_sim./sum(hist_sim);
            
            hold on;
            if exist('fitHistPlot','var') == 1
                delete(fitHistPlot)
            end
            fitHistPlot = plot(FRET_bins,hist_sim);
            fitHistPlot.LineStyle = '-';
            fitHistPlot.Color = 'red';
            fitHistPlot.LineWidth = 2;
            fitHistPlot.DisplayName = ['Gen' num2str(genNum)];
            
            if exist('state1_HistPlot','var') == 1
                delete(state1_HistPlot);
                delete(state2_HistPlot);
                delete(state3_HistPlot);
            end
            
            state1_HistPlot = plot(FRET_bins,p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2)./denom_hist_sim,'c--','LineWidth',1,'DisplayName','p1_{eq}');
            state2_HistPlot = plot(FRET_bins,p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2)./denom_hist_sim,'m--','LineWidth',1,'DisplayName','p2_{eq}');
            state3_HistPlot = plot(FRET_bins,p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2)./denom_hist_sim,'g--','LineWidth',1,'DisplayName','p3_{eq}');
            legend('show');
        end
        
        %--------------------------------------------------------------------------
        % (2) Optimization Target #2: 2-point TCF (C2)           (20 uSec
        % and on)  (PART 3: Make the % initial generation of guesses)
        %--------------------------------------------------------------------------
        if fitC2Mode == 1
            figure(2)
            
            hold on;
            if exist('C2plot','var') == 1
                delete(C2plot)
            end
            
            if useOldCodeMode == 1
                
                k12 = 1/t12;
                k13 = 1/t13;
                k21 = 1/t21;
                k23 = 1/t23;
                k31 = 1/t31;
                k32 = k12*k23*k31/(k13*k21);
                %         function tcf = TCF_cyclic3state(time,A0,A1,A2,k01,k10,k12,k21,k20)
                C2_sim = TCF_cyclic3state(C2_exp_x,A1,A2,A3,k12,k21,k23,k32,k31);
            elseif useOldCodeMode == 0
                C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
            end
            
            if normalizeMode == 1
                C2_sim = C2_sim./C2_sim(1);
            end
            C2plot = plot(C2_exp_x,C2_sim);
            set(gca,'xscale','log');
            drawnow();
        end
        
    end
    
    %--------------------------------------------------------------------------
    % Asses the chisquared of the guess
    %--------------------------------------------------------------------------
    tic
    if useOldCodeMode == 1
        k12 = 1/t12;
        k13 = 1/t13;
        k21 = 1/t21;
        k23 = 1/t23;
        k31 = 1/t31;
        
        chisquared = multigoaltcf_analytical_3state(k12,k13,k21,k23,k31,A1,A2,A3);
    else
        chisquared = multigoaltcf_analytical_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
    end
    pop_chisquared(pop_idx) = chisquared;
    
    %----------------------------------------------------------------------
    % (4) Plot the chisquared value of each member of the population
    %----------------------------------------------------------------------
    if plot_Gen1ChiSquaredMode == 1
        figure(4);
        
        set(gcf,'Name','Fitting Progress');
        if useFigurePosnMode == 1
            switch computer_terminal_str
                case 'computer_baimi_mode'
                    set(gcf,'Position',[1963 48 560 420]);%Work Desktop
                case 'computer_bisraels_mode'
                    set(gcf,'Position',[560 535 560 420]);%Macbook pro
            end
        end
        set(gcf,'Color','w');
        set(gca,'yscale','log');
        plot(pop_idx,chisquared,'--gs',...
            'LineWidth',2,...
            'MarkerSize',10,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[0.5,0.5,0.5])
        hold on;
        axis tight;
        xlim([1,inf]);
        
        xlabel('Generation One','FontSize',14);
        ylabel('RMS Fitness (\chi^2)','FontSize',14);
        title('\chi^2 of Generation 1','FontSize',14)
        drawnow();
    end
end


elapsedTime = toc;
if clockMode == 1
    disp(['     It took ' num2str(elapsedTime) ' seconds to compute the '...
        'RMS of ' num2str(NmembersInitPop) ' guesses.(' num2str(elapsedTime/NmembersInitPop) ' sec/guess)']);
end



%--------------------------------------------------------------------------
% (7) Make a histogram of the chisquared values
%--------------------------------------------------------------------------
if plot_ChiSquaredHistogramMode == 1
    figure(7)
    set(gcf,'Name','Chi-squared values of gen1 guesses');
    set(gcf,'Color','w');
    switch computer_terminal_str
        case 'computer_baimi_mode'
%             set(gcf,'Position',[2535 573 560 420]);%Work computer
         set(gcf,'Position',[1443 72 459 332]);%small widow
        case 'computer_bisraels_mode'
            %                 set(gcf,'Position',[1121 41 560 420]);%Macbook computer
            disp('Set an option for placement of this figure');
    end
    initialHist = histogram(pop_chisquared);
    xlabel('\chi^2 Value');
    ylabel('Frequency');
    title('Distribution of Gen1 guesses chi-squared');
end


%--------------------------------------------------------------------------
% Sort Generation 1 as a function of chisquare
%--------------------------------------------------------------------------
% Sort from Lowest rms to the highest rms (ascending order)\

[pop_chisquare_sorted,index_arr] = sort(pop_chisquared);%First in array have the lowest chisquare value


population = population(index_arr,:);

guess(genNum,1:8) = population(1,:); % Current best guess

Best_chisquared = pop_chisquare_sorted(1);
guess(genNum,9) = Best_chisquared; % Current best guess' fit value

genNum_array = genNum;
chisquared_array = Best_chisquared;

% Display best guess to screen

t12 = guess(genNum,1);
t13 = guess(genNum,2);
t21 = guess(genNum,3);
t23 = guess(genNum,4);
t31 = guess(genNum,5);

A1 = guess(genNum,6);
A2 = guess(genNum,7);
A3 = guess(genNum,8);

k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;

k32 = k12*k23*k31/(k13*k21);
t32 = 1/k32;

chisquared = guess(genNum,9);

% if verboseMode == 1
fprintf(['Best fit from initial generation was member #%d:\n t12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
    '\n A1 = %f, A2 = %f, A3 = %f, chisquared = %f\r\n'],...
    index_arr(1),t12,t13,t21,t23,t31,t32,A1,A2,A3,chisquared);
% end

%--------------------------------------------------------------------------
% (1) Make a histogram which corresponds to the best guess
%--------------------------------------------------------------------------
if fitHistMode == 1
    figure(1)
    hold on;
    
    if useOldCodeMode == 1
        k12 = 1/t12;
        k13 = 1/t13;
        k21 = 1/t21;
        k23 = 1/t23;
        k31 = 1/t31;
        k32 = k12*k23*k31/(k13*k21);
        % %         function [p0_eq,p1_eq,p2_eq,lam1,lam2] = FourPtTCF_cyclic3state_xo(tau1range,tau2,A0,A1,A2,k01,k10,k12,k21,k20)
        %         [p1_eq,p2_eq,p3_eq,~,~] = FourPtTCF_cyclic3state_xo(C4_tau1range,0,A1,A2,A3,k12,k21,k23,k32,k31);
        [p1_eq,p2_eq,p3_eq] = cyclic3state_hist(A1,A2,A3,k12,k21,k23,k32,k31);
    elseif useOldCodeMode == 0
        [Peq] = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
        p1_eq = Peq(1);
        p2_eq = Peq(2);
        p3_eq = Peq(3);
    end
    hist_sim = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim = sum(hist_sim);
    hist_sim = hist_sim./denom_hist_sim;
    
    hold on;
    if exist('fitHistPlot','var') == 1
        delete(fitHistPlot)
    end
    fitHistPlot = plot(FRET_bins,hist_sim);
    fitHistPlot.LineStyle = '-';
    fitHistPlot.Color = 'red';
    fitHistPlot.LineWidth = 2;
    fitHistPlot.DisplayName = ['Gen' num2str(genNum)];
    
    if exist('state1_HistPlot','var') == 1
        delete(state1_HistPlot);
        delete(state2_HistPlot);
        delete(state3_HistPlot);
    end
    
    state1_HistPlot = plot(FRET_bins,p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2)./denom_hist_sim,'c--','LineWidth',1,'DisplayName','p1_{eq}');
    state2_HistPlot = plot(FRET_bins,p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2)./denom_hist_sim,'m--','LineWidth',1,'DisplayName','p2_{eq}');
    state3_HistPlot = plot(FRET_bins,p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2)./denom_hist_sim,'g--','LineWidth',1,'DisplayName','p3_{eq}');
    legend('show');
end

%--------------------------------------------------------------------------
% (4) Make a red circle around the best fit and begin optimization from there
%--------------------------------------------------------------------------
if plot_Gen1ChiSquaredMode == 1
    figure(4);
    plot( index_arr(1),chisquared,'--go','LineWidth',3,...
        'MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','none')
    y = [chisquared chisquared];
    x = [index_arr(1) 1];
    line1 = line(x,y);
    line1.LineWidth = 2;
    line1.Color = 'r';
    line1.LineStyle = '--';
end


%--------------------------------------------------------------------------
% (4) Keep track of the chisquared to make sure it drops with each generation
%--------------------------------------------------------------------------
if plot_GenerationalProgressMode == 1
    figure(4);
    %     clf;
    set(gcf,'Name','Fitting Progress');
    if useFigurePosnMode == 1
        switch computer_terminal_str
            case 'computer_baimi_mode'
                set(gcf,'Position',[1963 48 560 420]);%Work Desktop
            case 'computer_bisraels_mode'
                set(gcf,'Position',[560 535 560 420]);%Macbook pro
        end
    end
    set(gcf,'Color','w');
    
    plot(genNum,chisquared,'bs',...
        'MarkerSize',10,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0.5,0.5,0.5]);
    
    axis tight;
    xlim([1,inf]);
    drawnow();
    
    xlabel('Member Number','FontSize',14);
    ylabel('RMS Fitness (\chi^2)','FontSize',14);
    title('Fitness vs Generation Number','FontSize',14)
    
end
%--------------------------------------------------------------------------
% (5) Plot histograms of initial guesses for each paramater
%--------------------------------------------------------------------------
figure(5)
set(gcf,'Name','Distribution of initial paramaters');
set(gcf,'Color','w');
switch computer_terminal_str
    case 'computer_baimi_mode'
        set(gcf,'Position',[2536 45 1260 425]);%Work computer
    case 'computer_bisraels_mode'
        set(gcf,'Position',[9 46 1111 425]);%Macbook computer
end
param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
maxCountArray = zeros(1,Nparam);
for param_idx = 1:8
    subplot(2,4,param_idx)
    
    if ismember(param_idx,1:5)%Rate histograms
        nbins = 50;
        histogram(population(:,param_idx),nbins);
        [N,edges] = histcounts(population(:,param_idx),nbins);
        maxCountArray(param_idx) = max(N);
        xlabel('Time (Sec)');
    elseif ismember(param_idx,[6,7,8]) == 1%FRET Historgrams
        nbins = 20;
        histogram(population(:,param_idx),nbins);
        [N,edges] = histcounts(population(:,param_idx),nbins);
        maxCountArray(param_idx) = max(N);
        xlabel('FRET');
        xlim([0,1]);
    end
    %Plot a blue line for the best guess of that gneration
    progressLine = line([guess(genNum,param_idx) guess(genNum,param_idx) ],[0  maxCountArray(param_idx)],...
        'Color','blue','LineStyle','-','LineWidth',4);
    
    ylabel('Frequency');
    title(cell2str(param_strings(param_idx)));
end

% error('Quitting before generation 2');
if pauseBetweenGenerationMode == 1
    disp(['     Gen' num2str(genNum) ' Chi^2 = ' num2str(chisquared) ': Press Enter to proceed to the next generation']);
    pause();
end


%//////////////////////////////////////////////////////////////////////////
% PART 4: Refine the paramaters by mutating certain paramaters randomly
%//////////////////////////////////////////////////////////////////////////
%Bookmarking for keeping up with the number of trials
Ntrials = 0;
%
for genNum = 2:maxGenerations
    tic
    % Mix genes of top rated individuals (population is sorted by chisquared: The lowest chisquared are first)
    for n = 1:NmembersInitPop - membersToKeep
        newpopulation(n,:) = [population(ceil(Nreproduce*rand(1,1)),1),population(ceil(Nreproduce*rand(1,1)),2),population(ceil(Nreproduce*rand(1,1)),3),population(ceil(Nreproduce*rand(1,1)),4),population(ceil(Nreproduce*rand(1,1)),5),population(ceil(Nreproduce*rand(1,1)),6),population(ceil(Nreproduce*rand(1,1)),7),population(ceil(Nreproduce*rand(1,1)),8)];
    end
    
    % Keep top 5% best parents (These are the largest indices)
%     newpopulation(NmembersInitPop - 1,:) = population(2,:);
%     newpopulation(NmembersInitPop,:) = population(1,:);
    pop_idx = 0;
    while pop_idx <= membersToKeep
%         Let the highest member#s be the best guesses from the first
%         generation
        newpopulation(NmembersInitPop - pop_idx,:) = population(pop_idx+1,:);
    pop_idx = pop_idx + 1;
    end
    
    
    
    population = newpopulation;
    
    % Mutate some individuals
    for q = 1:Nmutations
        %Let each paramater be mutated by a different amount across entire
        %population
        rndArr = rand(Nparam,1);
        for param_idx = 1:Nparam
            % Move each value of the population somewhere inbetween (uniformly
            % distributed) up or down by its  maximum allowable mutatin value
            popMemberID = ceil((NmembersInitPop - 2)*rndArr(param_idx));
            newValue = population(popMemberID,param_idx) + (2*(rand(1,1) - .5))*maxMutationArray(param_idx);
            population(popMemberID,param_idx) = newValue;
            if guessUpdateMode == 1
                fprintf(' Mutating member #%d''s %s param to %f\r',popMemberID,cell2str(param_strings(param_idx)),newValue)
            end
            % If the population has exceeded the maximum allowable amount, set it to
            % that amount. If it is below the minumum amont, set it to the minimum
            %************ How often does this happen? This is definitely not what
            %we want to occur*********
            if population(popMemberID,param_idx) > boundsArray(param_idx,2)
                population(popMemberID,param_idx) = boundsArray(param_idx,2);
                if guessUpdateMode == 1
                    fprintf('     Setting param %d (%s) to the MAX value: %f.\n',param_idx,char(param_strings(param_idx)),boundsArray(param_idx,2));
                end
            elseif population(popMemberID,param_idx) < boundsArray(param_idx,1)
                population(popMemberID,param_idx) = boundsArray(param_idx,1);
                if guessUpdateMode == 1
                    fprintf('     Setting param %d (%s) to the MIN value: %f.\n',param_idx,char(param_strings(param_idx)),boundsArray(param_idx,1));
                    maxMutationCountsArray(param_idx) = maxMutationCountsArray(param_idx) + 1;
                end
            else
                if guessUpdateMode == 1
                    fprintf('     Setting param %d (%s) to the NEW value: %f.\n',param_idx,char(param_strings(param_idx)),newValue);
                    minMutationCountsArray(param_idx) = minMutationCountsArray(param_idx) + 1;
                end
            end
        end
        
    end
    
    %----------------------------------------------------------------------
    %After Mutations, Recalculate the chisquared values of each member
    %----------------------------------------------------------------------
    for pop_idx = 1:NmembersInitPop
        %     parfor pop_idx = 1:NmembersInitPop
        %         popfit(k) = rmscalc(1/population(k,1),1/population(k,2),1/population(k,3),1/population(k,4),1/population(k,5),population(k,6),population(k,7),population(k,8));
        %     chisquared = multigoaltcf_analytical_3state(1/population(k,1),1/population(k,2),1/population(k,3),1/population(k,4),1/population(k,5),population(k,6),population(k,7),population(k,8));
        
        
        t12 = population(pop_idx,1);
        t13 = population(pop_idx,2);
        t21 = population(pop_idx,3);
        t23 = population(pop_idx,4);
        t31 = population(pop_idx,5);
        
        A1 = population(pop_idx,6);
        A2 = population(pop_idx,7);
        A3 = population(pop_idx,8);
        
        
        if useOldCodeMode == 1
            k12 = 1/t12;
            k13 = 1/t13;
            k21 = 1/t21;
            k23 = 1/t23;
            k31 = 1/t31;
            chisquared = multigoaltcf_analytical_3state(k12,k13,k21,k23,k31,A1,A2,A3);
        else
            chisquared = multigoaltcf_analytical_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
        end
        %Add the chisquared value to its corresponding slot in the
        %population
        pop_chisquared(pop_idx) = chisquared;
    end
    
    % Evaluate Generation
    % Sort from highest rated to lowest rated individuals
    [pop_chisquare_sorted,index_arr] = sort(pop_chisquared);
    population = population(index_arr,:);
    
    %Add to the best guess array the best guess of the population
    guess(genNum,1:8) = population(1,:); % Current best guess
    guess(genNum,9) = pop_chisquare_sorted(1); % Current best guess' fit value
    
    %population(1,:) % Display best guess to screen
    %     popfit(1);
    if verboseMode == 1
        disp(['     Gen ' num2str(genNum) ': Best fitness guess so far: ' num2str(pop_chisquare_sorted(1))]);
    end
    %--------------------------------------------------------------------------
    % (6) Display a graphic of the network (PART 4: genX)
    %--------------------------------------------------------------------------
    figure(6);
    clf
    switch computer_terminal_str
        case 'computer_baimi_mode'
            set(gcf,'Position',[3113 571 560 420]);%Work Desktop
        case 'computer_bisraels_mode'
            set(gcf,'Position',[1121 535 560 420]);%Macbook pro
    end
    set(gcf,'Name',['Gen' num2str(genNum) ': Model Results']);
    
    set(gcf,'Color',[1 1 1]);
    hold on;
    xlim([0 1.1]);
    ylim([0 1]);
    
    %Designate spots for the states
    state1_loc = [0.5 1];
    state2_loc = [1 0];
    state3_loc = [0 0];
    
    %Plot the gen number
    text(state1_loc(1) + (state2_loc(1) - state1_loc(1))/2,state1_loc(2),['Gen# = ' num2str(genNum)],'FontSize',14);
    
    %Plot the state symbols
    text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
    text(state2_loc(1),state2_loc(2),'2','FontSize',24);
    text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
    
    %Plot the FRET values
    text(state1_loc(1),state1_loc(2) + 0.05,['=' num2str(A1,'%.2f')],'FontSize',16)
    text(state2_loc(1)+ 0.05,state2_loc(2),['=' num2str(A2,'%.2f')],'FontSize',16);
    text(state3_loc(1),state3_loc(2),['=' num2str(A3,'%.2f')],'FontSize',16);
    
    %Plot Lines Between states
    line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
    line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
    line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
    
    %Plot the inverse of the rates
    t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
    t12_msg = ['t_{12} = ' num2str(t12) 10 't_{21} = ' num2str(t21)];
    text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',14);
    
    t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
    t23_msg = ['t_{23} = ' num2str(t23) 10 't_{32} = ' num2str(t32)];
    text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',14);
    
    t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
    t31_msg = ['t_{31} = ' num2str(t31) 10 't_{13} = ' num2str(t13)];
    text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',14);
    
    %Plot the chi-squared value
    text(state3_loc(1),state1_loc(2),['\chi^2 = ' num2str(chisquared)],'FontSize',14);
    
    
    axis off;
    
    %----------------------------------------------------------------------
    % (7) Plot the distribution of guesses (PART 4: genX)
    %----------------------------------------------------------------------
    if plot_ChiSquaredHistogramMode == 1
        figure(7)
        set(gcf,'Name','Chi-squared values of initial guesses');
        set(gcf,'Color','w');
%         switch computer_terminal_str
%             case 'computer_baimi_mode'
%                 set(gcf,'Position',[2535 573 560 420]);%Work computer
%             case 'computer_bisraels_mode'
%                 %                 set(gcf,'Position',[1121 41 560 420]);%Macbook computer
%                 disp('Set an option for placement of this figure');
%         end
        histogram(pop_chisquared);
        xlabel('\chi^2 Value');
        ylabel('Frequency');
        title(['Gen' num2str(genNum) ' Distribution of guesses chi-squared']);
    end
    
    %----------------------------------------------------------------------
    % Determine if the guesses have converged (PART 4: genX)
    %----------------------------------------------------------------------
    if (guess(genNum-1,9) - guess(genNum,9))/guess(genNum,9) <= threshold
        Ntrials = Ntrials + 1;
    else
        Ntrials = 0;
    end
    if Ntrials == maxRepeats
        guess = guess(1:genNum,:);
        disp(['The fitness has no longer improved after ' num2str(maxRepeats) ' repeats. Ending search.']);
        break;
    end
    
    if clockMode == 1
        elapsedTime = toc;
        disp(['The amount of time per generation is: ' num2str(elapsedTime) ' seconds']);
    end
    
    %--------------------------------------------------------------------------
    % Update the plots with the current best guesses (PART 4: genX)
    %--------------------------------------------------------------------------
    if plot_GenerationalProgressMode == 1
        
        t12 = guess(genNum,1);
        t13 = guess(genNum,2);
        t21 = guess(genNum,3);
        t23 = guess(genNum,4);
        t31 = guess(genNum,5);
        
        A1 = guess(genNum,6);
        A2 = guess(genNum,7);
        A3 = guess(genNum,8);
        
        chisquared = guess(genNum,9);
        
        genNum_array = [genNum_array; genNum];
        chisquared_array = [chisquared_array; chisquared];
        
        %--------------------------------------------------------------------------
        % (1) Optimization Target #1: 1-D FRET HISTOGRAM (PART 4: genX)
        %--------------------------------------------------------------------------
        if fitHistMode == 1
            figure(1)
            hold on;
            
            if useOldCodeMode == 1
                
                k12 = 1/t12;
                k13 = 1/t13;
                k21 = 1/t21;
                k23 = 1/t23;
                k31 = 1/t31;
                k32 = k12*k23*k31/(k13*k21);
                % %         function [p0_eq,p1_eq,p2_eq,lam1,lam2] = FourPtTCF_cyclic3state_xo(tau1range,tau2,A0,A1,A2,k01,k10,k12,k21,k20)
                %                 [p1_eq,p2_eq,p3_eq,~,~] = FourPtTCF_cyclic3state_xo(C4_tau1range,0,A1,A2,A3,k12,k21,k23,k32,k31);
                [p1_eq,p2_eq,p3_eq] = cyclic3state_hist(A1,A2,A3,k12,k21,k23,k32,k31);
            else
                [Peq] = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
                p1_eq = Peq(1);
                p2_eq = Peq(2);
                p3_eq = Peq(3);
            end
            
            hist_sim = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
            denom_hist_sim = sum(hist_sim);
            hist_sim = hist_sim./sum(hist_sim);
            
            hold on;
            if exist('fitHistPlot','var') == 1
                delete(fitHistPlot)
            end
            fitHistPlot = plot(FRET_bins,hist_sim);
            fitHistPlot.LineStyle = '-';
            fitHistPlot.Color = 'red';
            fitHistPlot.LineWidth = 2;
            fitHistPlot.DisplayName = ['Gen' num2str(genNum)];
            
            
            if exist('state1_HistPlot','var') == 1
                delete(state1_HistPlot);
                delete(state2_HistPlot);
                delete(state3_HistPlot);
            end
            
            state1_HistPlot = plot(FRET_bins,p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2)./denom_hist_sim,'c--','LineWidth',1,'DisplayName','p1_{eq}');
            state2_HistPlot = plot(FRET_bins,p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2)./denom_hist_sim,'m--','LineWidth',1,'DisplayName','p2_{eq}');
            state3_HistPlot = plot(FRET_bins,p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2)./denom_hist_sim,'g--','LineWidth',1,'DisplayName','p3_{eq}');
            legend('show');
            drawnow();
            
        end
        
        %--------------------------------------------------------------------------
        % (2) Optimization Target #2: 2-point TCF (C2) (20 uSec and on) (PART 4: genX)
        %--------------------------------------------------------------------------
        if fitC2Mode == 1
            figure(2)
            
            hold on;
            if exist('C2plot','var') == 1
                delete(C2plot)
            end
            
            if useOldCodeMode == 1
                
                k12 = 1/t12;
                k13 = 1/t13;
                k21 = 1/t21;
                k23 = 1/t23;
                k31 = 1/t31;
                k32 = k12*k23*k31/(k13*k21);
                %         function tcf = TCF_cyclic3state(time,A0,A1,A2,k01,k10,k12,k21,k20)
                C2_sim = TCF_cyclic3state(C2_exp_x,A1,A2,A3,k12,k21,k23,k32,k31);
            else
                C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
            end
            
            if normalizeMode == 1
                C2_sim = C2_sim./C2_sim(1);
            end
            C2plot = plot(C2_exp_x,C2_sim);
            set(gca,'xscale','log');
            drawnow();
        end
        
        %--------------------------------------------------------------------------
        % (3) Optimization Target #3: 4-point TCF (C4) (10 uSec and on) (PART 4: genX)
        %--------------------------------------------------------------------------
        if fitC4Mode == 1
        end
        
        %--------------------------------------------------------------------------
        % (4) Plotting the progress as a function of Generation
        %--------------------------------------------------------------------------
        figure(4);
        if useFigurePosnMode == 1
            switch computer_terminal_str
                case 'computer_baimi_mode'
                    set(gcf,'Position',[1963 48 560 420]);%Work Desktop
                case 'computer_bisraels_mode'
                    set(gcf,'Position',[560 535 560 420]);%Macbook pro
            end
        end
        
        hold on;
        plot(genNum_array,chisquared_array,'o',...
            'MarkerSize',10,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        drawnow();
        %--------------------------------------------------------------------------
        % (5) Update paramter-histograms as a function of the current best
        % guesses  (PART 4: genX)
        %--------------------------------------------------------------------------
        figure(5)
        set(gcf,'Name','Distribution of initial paramaters');
        set(gcf,'Color','w');
        switch computer_terminal_str
            case 'computer_baimi_mode'
                set(gcf,'Position',[2536 45 1260 425]);%Work computer
            case 'computer_bisraels_mode'
                set(gcf,'Position',[9 46 1111 425]);%Macbook computer
        end
        param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
        maxCountArray = zeros(1,Nparam);
        for param_idx = 1:8
            subplot(2,4,param_idx)
            
            if ismember(param_idx,1:5)%Rate histograms
                nbins = 50;
                histogram(population(:,param_idx),nbins);
                [N,edges] = histcounts(population(:,param_idx),nbins);
                maxCountArray(param_idx) = max(N);
                xlabel('Time (Sec)');
            elseif ismember(param_idx,[6,7,8]) == 1%FRET Historgrams
                nbins = 20;
                histogram(population(:,param_idx),nbins);
                [N,edges] = histcounts(population(:,param_idx),nbins);
                maxCountArray(param_idx) = max(N);
                xlabel('FRET');
                xlim([0,1]);
            end
            if exist('progressLine','var') == 1
                delete(progressLine)
            end
            %Plot a blue line for the best guess of that gneration
            progressLine = line([guess(genNum,param_idx) guess(genNum,param_idx) ],[0  maxCountArray(param_idx)],...
                'Color','blue','LineStyle','-','LineWidth',4);
            
            ylabel('Frequency');
            title(cell2str(param_strings(param_idx)));
        end
        
    end
    
    if pauseBetweenGenerationMode == 1
        disp(['     Gen' num2str(genNum) ' Chi^2 = ' num2str(chisquared) ': Press Enter to proceed to the next generation']);
        pause();
    end
end

%//////////////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------
% PART 5: Display the Final values to the user
%--------------------------------------------------------------------------
%//////////////////////////////////////////////////////////////////////////

t12 = guess(genNum,1);
t13 = guess(genNum,2);
t21 = guess(genNum,3);
t23 = guess(genNum,4);
t31 = guess(genNum,5);

A1 = guess(genNum,6);
A2 = guess(genNum,7);
A3 = guess(genNum,8);

k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;

k32 = k12*k23*k31/(k13*k21);
t32 = 1/k32;

chisquared = guess(genNum,9);
newChisquared = chisquared;
fprintf(['Final values before checking old fits:\n t12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
    '\n A1 = %f, A2 = %f, A3 = %f, chisquared = %f\r\n'],...
    t12,t13,t21,t23,t31,t32,A1,A2,A3,chisquared);


%--------------------------------------------------------------------------
% Save the data (PART 5: Final Values)
%--------------------------------------------------------------------------
if saveMode == 1
    %--------------------------------------------------------------------------
    % Make a folder to hold output (PART 5: Final Values)
    %--------------------------------------------------------------------------
    
    %Make output folder if it doesnt exist
    outputFolderName = [programName '_output'];
    if exist(outputFolderName,'dir') ~= 7
        mkdir(outputFolderName);
        disp('Making a folder to hold the output');
    end
    
    if exist([outputFolderName filesep() 'BestFitResults' '.mat'],'file') == 2
        load([outputFolderName filesep() 'BestFitResults' '.mat'],'chisquared','iter');
        oldChisquared = chisquared;
        if newChisquared < oldChisquared
            iter = iter+1;
            fprintf('Congratulations! %f < %f. Updating best fit.\r\n',newChisquared,oldChisquared)
            chisquared = newChisquared;
            %Save the results of the fitting routine
            
            disp(['Saving the filenames in ' outputFolderName filesep() 'fitInputData.mat']);
            
            if fitHistMode == 1
                save([outputFolderName filesep() 'fitInputData.mat'],'histogram_FileName','FRET_bins','targetHistogram','-append')
            end
            if fitC2Mode == 1
                save([outputFolderName filesep() 'fitInputData.mat'],'C2_FileName','C2_exp_x','C2_exp_y','weightC2func','-append');
            end
            if fitC4Mode == 1
                save([outputFolderName filesep() 'fitInputData.mat'],'C4_FileName','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','-append');
            end
            
            save([outputFolderName filesep() 'BestFitResults' '.mat'],'A1','A2','A3',...
                't12','t13','t21','t23','t31','t32','chisquared','guess',...
                'k12','k13','k21','k23','k31','k32','iter');
            
            %Save the paramaters used for the genetic algorithm
            save([outputFolderName filesep() 'genAlgParamaters.mat'],'programName','threshold',...
                'NmembersInitPop','maxGenerations','maxRepeats','Nreproduce',...
                'Nmutations','boundsArray','maxMutationArray')
            
            %Save paramatres used to plot
            save([outputFolderName filesep() 'plottingParamaters.mat'],'FRET_bins','sigma_A1','sigma_A2','sigma_A3',...
                'C2_exp_x','C2_exp_y');
            
        elseif newChisquared == oldChisquared
            iter = iter + 1;
            save([outputFolderName filesep() 'BestFitResults'],'iter','-append');
            fprintf('Nothing changed.%f = %f.  Updated the BestFitResults with another iteration. \r\n',newChisquared,oldChisquared);
            
        elseif newChisquared > oldChisquared
            saveMode = 0; %no need to save anything going forward
            iter = iter + 1;
            save([outputFolderName filesep() 'BestFitResults'],'iter','-append');
            fprintf('Oh No! %f > %f. Keeping old fit. Updated the BestFitResults with another iteration. \r\n',newChisquared,oldChisquared)
            if showBestFitMode == 1
                disp('     Loading the best fit from a previous run.');
                load([outputFolderName filesep() 'BestFitResults.mat'],'A1','A2','A3',...
                    't12','t13','t21','t23','t31','t32','chisquared','guess',...
                    'k12','k13','k21','k23','k31','k32','iter');
            end
        end
    else
        %This is the case where you have found no best fit file yet.
        fprintf('Saving the fit for the first time. chisquare = %f\r',newChisquared);
        iter = 0;
        %Save the results of the fitting routine
        save([outputFolderName filesep() 'BestFitResults' '.mat'],'A1','A2','A3',...
            't12','t13','t21','t23','t31','t32','chisquared','guess',...
            'k12','k13','k21','k23','k31','k32','iter');
        
        %Save the paramaters used for the genetic algorithm
        save([outputFolderName filesep() 'genAlgParamaters.mat'],'programName','threshold',...
            'NmembersInitPop','maxGenerations','maxRepeats','Nreproduce',...
            'Nmutations','boundsArray','maxMutationArray')
        
        %Save paramatres used to plot
        save([outputFolderName filesep() 'plottingParamaters.mat'],'FRET_bins','sigma_A1','sigma_A2','sigma_A3',...
            'C2_exp_x','C2_exp_y');
        
    end
    fprintf(['Final values selected:\n t12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
        '\n A1 = %f, A2 = %f, A3 = %f, chisquared = %f\r\n'],...
        t12,t13,t21,t23,t31,t32,A1,A2,A3,chisquared);
    
end



%//////////////////////////////////////////////////////////////////////////
% Plot final best-fit results of the algorithm
%//////////////////////////////////////////////////////////////////////////

%--------------------------------------------------------------------------
% (1) Optimization Target #1: 1-D FRET HISTOGRAM (PART 5: Final Values)
%--------------------------------------------------------------------------
if plotMode == 1
    figure(1)
    
    set(gcf,'Name','FRET Histogram');
    if useFigurePosnMode == 1
        switch computer_terminal_str
            case 'computer_baimi_mode'
                set(gcf,'Position',[1963 573 560 420]);%Work Desktop
            case 'computer_bisraels_mode'
                set(gcf,'Position',[560 535 560 420]);%Macbook pro
        end
        
    end
    set(gcf,'Color','w');
    %Calculate the equilibrium Populations
    if useOldCodeMode == 0
        [Peq] = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
        p1_eq = Peq(1);
        p2_eq = Peq(2);
        p3_eq = Peq(3);
    elseif useOldCodeMode == 1
        
        %     function [p0_eq,p1_eq,p2_eq,lam1,lam2] = FourPtTCF_cyclic3state_xo(tau1range,tau2,A0,A1,A2,k01,k10,k12,k21,k20)
        [p1_eq, p2_eq, p3_eq, lam1, lam2] = FourPtTCF_cyclic3state_xo(1:2,0,A1,A2,A3,k12,k21,k23,k32,k31);
        disp(['The final eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/lam1) ' microseconds and tau2 = 1/lam2 = ' num2str(1e6*1/lam2) ' microseconds']);
    end
    
    hist_sim_final = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim_final = sum(hist_sim_final);
    %     if exist('normalizeMode','var') == 1
    % %         if normalizeMode == 1
    %             hist_sim_final = hist_sim_final./denom_hist_sim_final;
    % %         end
    %     end
    hist_sim_final = hist_sim_final./denom_hist_sim_final;
    
    hold on;
    if exist('fitHistPlot','var') == 1
        delete(fitHistPlot)
    end
    
    fitHistPlot = plot(FRET_bins,hist_sim_final,'r-','LineWidth',2,'DisplayName','Final Fit (old code)');
    
    
    %---------------------------------------------------------------------
    % Plot the underlying FRET States (PART 5: Final Values)
    %---------------------------------------------------------------------
    if exist('state1_HistPlot','var') == 1
        delete(state1_HistPlot);
        delete(state2_HistPlot);
        delete(state3_HistPlot);
    end
    
    state1_HistPlot = plot(FRET_bins,p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2)./denom_hist_sim_final,'c--','LineWidth',1,'DisplayName','p1_{eq}');
    state2_HistPlot = plot(FRET_bins,p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2)./denom_hist_sim_final,'m--','LineWidth',1,'DisplayName','p2_{eq}');
    state3_HistPlot = plot(FRET_bins,p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2)./denom_hist_sim_final,'g--','LineWidth',1,'DisplayName','p3_{eq}');
    legend('show');
    
    
    if saveMode == 1
        saveas(gcf,[outputFolderName filesep()  'BestFitResults_hist.fig'])
    end
    
    %--------------------------------------------------------------------------
    % (2) Optimization Target #2: 2-point TCF (C2) (20 uSec and on) (PART 5: Final Values)
    %--------------------------------------------------------------------------
    if fitC2Mode == 1
        if useOldCodeMode == 1
            C2_sim = TCF_cyclic3state(C2_exp_x,A1,A2,A3,k12,k21,k23,k32,k31);
            display_str = 'TCF_cyclic3state';
        else
            C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
            display_str = 'C2Maker_3state123_cyclical';
        end
        figure(2)
        set(gcf,'Name','Final fit of data');
        if normalizeMode == 1
            C2_sim = C2_sim./C2_sim(1);
        end
        
        
        hold on;
        if exist('C2plot','var') == 1
            delete(C2plot)
        end
        C2plot = plot(C2_exp_x,C2_sim,'r-','LineWidth',2,'DisplayName',display_str);
    
    end
    
    
    %--------------------------------------------------------------------------
    % (3) Optimization Target #3: 4-point TCF (C4) (10 uSec and on) (PART 5: Final Values)
    %--------------------------------------------------------------------------
    if fitC4Mode == 1
        disp('Need code to display final 4point!!!!')
    end
    %--------------------------------------------------------------------------
    % (6) Display a graphic of the network (PART 5: Final Values)
    %--------------------------------------------------------------------------
    figure(6);
    clf
    switch computer_terminal_str
        case 'computer_baimi_mode'
            set(gcf,'Position',[3113 571 560 420]);%Work Desktop
        case 'computer_bisraels_mode'
            set(gcf,'Position',[1121 535 560 420]);%Macbook pro
    end
    
    set(gcf,'Name','Model Results');
    
    set(gcf,'Color',[1 1 1]);
    hold on;
    xlim([0 1.1]);
    ylim([0 1]);
    
    %Designate spots for the states
    state1_loc = [0.5 1];
    state2_loc = [1 0];
    state3_loc = [0 0];
    
    %Plot the state symbols
    text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
    text(state2_loc(1),state2_loc(2),'2','FontSize',24);
    text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
    
    %Plot the FRET values
    text(state1_loc(1),state1_loc(2) + 0.05,['=' num2str(A1,'%.2f')],'FontSize',16)
    text(state2_loc(1)+ 0.05,state2_loc(2),['=' num2str(A2,'%.2f')],'FontSize',16);
    text(state3_loc(1),state3_loc(2),['=' num2str(A3,'%.2f')],'FontSize',16);
    
    %Plot Lines Between states
    line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
    line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
    line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
    
    %Plot the inverse of the rates
    t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
    t12_msg = ['t_{12} = ' num2str(t12) 10 't_{21} = ' num2str(t21)];
    text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',14);
    
    t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
    t23_msg = ['t_{23} = ' num2str(t23) 10 't_{32} = ' num2str(t32)];
    text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',14);
    
    t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
    t31_msg = ['t_{31} = ' num2str(t31) 10 't_{13} = ' num2str(t13)];
    text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',14);
    
    %Plot the chi-squared value
    text(state3_loc(1),state1_loc(2),['\chi^2 = ' num2str(chisquared)],'FontSize',14);
    
    
    axis off;
    if saveMode == 1
        saveas(gcf,[outputFolderName filesep() 'ModelResultsFigure.fig'])
    end
end




%//////////////////////////////////////////////////////////////////////////
% OPTIMIZATION FUNCTION: multigoaltcf_analytical_3state
%//////////////////////////////////////////////////////////////////////////
% FUNCTION: multigoaltcf_analytical_3state
% PURPOSE: returns the rms (root mean square) deviation from the model
% INPUT: (1) FRET Values for each state
%        (2) Set of rates {k_ij}
function chisquared = multigoaltcf_analytical_3state(k12,k13,k21,k23,k31,A1,A2,A3)

global verboseMode normalizeMode diagnoseMode
global genNum targetHistogram weightingFactor_FREThist FRET_bins
global C2_exp_x C2_exp_y weightingFactor_C2 weightC2func
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
global fitHistMode fitC2Mode fitC4Mode
global sigma_A1 sigma_A2 sigma_A3
global yoff
chisquared_array = zeros(3,1);

%Calculate the rate given by detailed balance
k32 = k12*k23*k31/(k13*k21);
%------------------------------------------------------------------
% (1) Optimization Target #1: Get a single value for the entire
% hitogram: rms_array(1) (OPTIMIZATION FUNCTION:% multigoaltcf_analytical_3state)
%------------------------------------------------------------------
if fitHistMode == 1
 [p1_eq,p2_eq,p3_eq] = cyclic3state_hist(A1,A2,A3,k12,k21,k23,k32,k31);
    hist_sim = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
    hist_sim = hist_sim/sum(hist_sim);
  chisquared_array(1) = sum((hist_sim-targetHistogram ).^2)*weightingFactor_FREThist;
end

%--------------------------------------------------------------------------
% (2) Optimization Target #2: Get a single value for the entire 2-point TCF: rms_array(1)
%--------------------------------------------------------------------------
if fitC2Mode == 1
    % Calculate 2-pt. TCF using analytical formula
    %     function tcf = TCF_cyclic3state(time,A0,A1,A2,k01,k10,k12,k21,k20)
    C2_sim = TCF_cyclic3state(C2_exp_x,A1,A2,A3,k12,k21,k23,k32,k31);
%     C2_sim = C2_sim + yoff;

    if normalizeMode == 1
        C2_sim = C2_sim./C2_sim(1);
    end
    % Calculate rms deviation between C2_experiment and C2_simulated
    chisquared_array(2) = sum(((C2_sim-C2_exp_y).*weightC2func).^2)*weightingFactor_C2;
      
end
%--------------------------------------------------------------------------
% Optimization Target #3: Get a single value for the entire 4-point TCF: rms_array(3)
%--------------------------------------------------------------------------
% Calculate 4-pt. TCF from analytical formula with various tau2 times.
if fitC4Mode == 1
    [C4_tau2eq0_sim,~] = FourPtTCF_cyclic3state_norm(C4_tau1range,0,A1,A2,A3,k12,k21,k23,k32,k31);
    if normalizeMode == 1
        C4_tau2eq0_sim = C4_tau2eq0_sim./C4_tau2eq0_sim(1,1);
    end
    % Calculate rms deviation from 4-pt. TCF of guess with t2 = 0
    chisquared_array(3) = mean(mean((C4_tau2eq0_sim-C4_tau2eq0_exp).^2.*wC4func))*weightingFactor_C4_t0;
end
%--------------------------------------------------------------------------
% Sum all of the indivual root-mean-squares to get an overall rms
%--------------------------------------------------------------------------
chisquared = sum(chisquared_array);
% if verboseMode
%     if genNum > 1
%         disp(['Gen ' num2str(genNum) ' : The total rms is ' num2str(rms)]);
%     end
% end

%------------------------------------------------------------------
% Display relative contributions
%------------------------------------------------------------------
if diagnoseMode == 1
    if genNum > 1
        if fitHistMode == 1
            disp(['     The contribution of the Histogram to the rms_array is ' num2str(chisquared_array(1)*100/chisquared) '%']);
        end
        if fitC2Mode == 1
            disp(['     The contribution of the 2ptTCF to the rms_array is ' num2str(chisquared_array(2)*100/chisquared) '%']);
        end
        if fitC4Mode == 1
            disp(['     The contribution of the 4-point TCF to the rms_array is ' num2str(chisquared_array(3)*100/chisquared) '%']);
        end
    end
end

end

%//////////////////////////////////////////////////////////////////////////
% OPTIMIZATION FUNCTION: MINIMIZE chisquared
%//////////////////////////////////////////////////////////////////////////
% FUNCTION: multigoaltcf_analytical_3state123_cyclical
% PURPOSE: returns the rms (root mean square) deviation from the model
% INPUT: (1) FRET Values for each state
%        (2) Set of rates {k_ij}
% function chisquared = multigoaltcf_analytical_3state(k12,k13,k21,k23,k31,A1,A2,A3)
function chisquared = multigoaltcf_analytical_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3)
global verboseMode normalizeMode diagnoseMode
global genNum targetHistogram weightingFactor_FREThist FRET_bins
global C2_exp_x C2_exp_y weightingFactor_C2 weightC2func
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
global fitHistMode fitC2Mode fitC4Mode
global sigma_A1 sigma_A2 sigma_A3

chisquared_array = zeros(3,1);

%Calculate the rate given by detailed balance
% k32 = k12*k23*k31/(k13*k21);
%------------------------------------------------------------------
% (1) Optimization Target #1: Get a single value for the entire
% hitogram: rms_array(1)
%------------------------------------------------------------------
if fitHistMode == 1
 Peq = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
    P1EQ = Peq(1);
    P2EQ = Peq(2);
    P3EQ = Peq(3);
    
    hist_sim = P1EQ*exp(-((FRET_bins-A1)/sigma_A1).^2) + P2EQ*exp(-((FRET_bins-A2)/sigma_A2).^2) + P3EQ*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim = sum(hist_sim);
    hist_sim = hist_sim./denom_hist_sim;
  
    chisquared_array(1) = sum((hist_sim-targetHistogram ).^2)*weightingFactor_FREThist;

end

%--------------------------------------------------------------------------
% (2) Optimization Target #2: Get a single value for the entire 2-point TCF: rms_array(1)
%--------------------------------------------------------------------------
if fitC2Mode == 1
  C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
    
    if normalizeMode == 1
        C2_sim = C2_sim./C2_sim(1);
    end
    % Calculate rms deviation between C2_experiment and C2_simulated
    chisquared_array(2) = sum(((C2_sim-C2_exp_y).*weightC2func).^2)*weightingFactor_C2;
    
end
%--------------------------------------------------------------------------
% Optimization Target #3: Get a single value for the entire 4-point TCF: rms_array(3)
%--------------------------------------------------------------------------
% Calculate 4-pt. TCF from analytical formula with various tau2 times.
if fitC4Mode == 1
    %     [C4_tau2eq0_sim,~] = FourPtTCF_cyclic3state_norm(C4_tau1range,0,A1,A2,A3,k12,k21,k23,k32,k31);
    if normalizeMode == 1
        C4_tau2eq0_sim = C4_tau2eq0_sim./C4_tau2eq0_sim(1,1);
    end
    % Calculate rms deviation from 4-pt. TCF of guess with t2 = 0
    chisquared_array(3) = mean(mean((C4_tau2eq0_sim-C4_tau2eq0_exp).^2.*wC4func))*weightingFactor_C4_t0;
    
end
%--------------------------------------------------------------------------
% Sum all of the indivual root-mean-squares to get an overall rms
%--------------------------------------------------------------------------
chisquared = sum(chisquared_array);

%------------------------------------------------------------------
% Display relative contributions
%------------------------------------------------------------------
if diagnoseMode == 1
    if genNum > 1
        if fitHistMode == 1
            disp(['     The contribution of the Histogram to the rms_array is ' num2str(chisquared_array(1)*100/chisquared) '%']);
        end
        if fitC2Mode == 1
            disp(['     The contribution of the 2ptTCF to the rms_array is ' num2str(chisquared_array(2)*100/chisquared) '%']);
        end
        if fitC4Mode == 1
            disp(['     The contribution of the 4-point TCF to the rms_array is ' num2str(chisquared_array(3)*100/chisquared) '%']);
            endc
        end
    end
    
end
end
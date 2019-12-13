%__________________________________________________________________________
% AUTHOR: Claire Albrecht and Bress Israels (November 2019)
%
% NAME: GenAlg_6state.m
%
% FUNCTION: % Fits (1) the FRET Histogram, (2) the 2pt TCF and (3) the 4-point TCF
%
%--------------------------------------------------------------------------
% PROCEDURE:
%--------------------------------------------------------------------------
% (1) Load the experimental histogram
% (2) Load the experimental 2pt time correlation functions to fit
% (3) Load the experimental 4pt TCFs for various tau2 values
% (4) Make an initial population of guesses for the model paramaters
% (5) Construct Histograms, C2, and {C4} for each guess
% (6) Compare each simulated plot to the experimental one to calculate rms
% (7) Choose the best elements of each and mix together.
% (8) Repeat steps 5-7 until the simulation converges below the threshold
%
%--------------------------------------------------------------------------
% INPUT:     !! NEED THESE !!
%--------------------------------------------------------------------------
% (1) FRET Histogram : 
% (2) Two-Point TCF : 
% (3) Four-Point TCF : 
% (4) Solutions to Diff Eqns: 
% OUTPUT:
% (1) BestFitResults.mat         %All the fitting paramaters
% (2) BestFitRestults_hist.fig   %
% (3) fitInputData.mat
% (4) genAlgParamaters.mat
% (5) ModelResultsFigure.fig
% (6) plottingParamaters.mat
%
%--------------------------------------------------------------------------
% EXTERNAL PROGRAMS CALLS: (needs these codes in MATLAB PATH to function)
%--------------------------------------------------------------------------
% ANALYTICAL Algorithms (fast)
% (1) histMaker_3state123_cyclical_analytical    % Calculates Peq
% (2) C2maker_3state123_cyclical_analytical      % Calculates C2
% (3) C4maker_3state123_cyclical_analytical      % Calculates C4

% NUMERICAL Algorithms (slow)
% (0) ODEsolver_3state123_cyclical.m    % Calculates the conditional probabilities (symCondProb_3state123_cyclical.mat)
% (1) histMaker_3state123_cyclical      % Calculates Peq
% (2) C2Maker_3state123_cyclical        % Calculates C2
% (3) C4Maker_3state123_cyclical        % Calculates C4
%
%__________________________________________________________________________
programName = 'GenAlg_8state';
close all

%//////////////////////////////////////////////////////////////////////////
% PART 1: Set the program up with various options
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
modelName = '8state';
protein_str = 'gp32_0p5uM';
% protein_str = '0p5uMgp32';

% constructFolderNames = {'S4S5'};
constructFolderNames = {'S1S2'};
% constructFolderNames = {'S18S2'};
% constructFolderNames =  {'S19S5'};
% constructFolderNames = {'S19S5';'S1S2';'S4S5';'S18S2';};
%--------------------------------------------------------------------------
% Genetic Algorithm Patamaters
%--------------------------------------------------------------------------
NmembersInitPop = 80;          % Size of inital pop
maxGenerations = 1000;
maxRepeats = 10;
percentReproduce = 25;          % top 25% combine together (can be mutated)
percentToKeep = 50;             % Number of population to pass to next gen (can be mutated)
Nmutations = 1*NmembersInitPop;  % Number of mutations per generation
threshold = 0.001;              % Minimal allowable percentage  difference before program quits.
maxIterations = 10; % 1006;           % How is this different than maxGenerations?
forceMoreIterationsMode = 1;    % Independent iterations
NumberToExtend = 1;             
% NOTES: lowering the percent that reproduces but raising the percent to
% keep as high seems to help.
    % Why? What are these really doing in the GenAlg process?
    
%-------------------DERIVED PARAMATERS------------------
Nreproduce = round(NmembersInitPop * percentReproduce/100);  % # individuals who reproduce per generation.
membersToKeep = NmembersInitPop * percentToKeep/100;

%--------------------------------------------------------------------------
% Declare global variables
%--------------------------------------------------------------------------
global normalizeMode verboseMode guessUpdateMode diagnoseMode plotMode %useAnalyticalAlgorithmsMode 
global genNum fitHistMode fitC2Mode fitC4Mode
global targetHistogram weightingFactor_FREThist FRET_bins sigma_A1 sigma_A2 sigma_A3 sigma_A4 sigma_A5 sigma_A6 sigma_A7 sigma_A8 sigma_A9  % Histogram optimization
global C2_exp_x C2_exp_y weightingFactor_C2  weightC2func yoff
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func zoff

global  K %sigma_A P

%--------------------------------------------------------------------------
% User Options                                                  (Part 1)
%--------------------------------------------------------------------------
normalizeMode = 0;      % Normalizes the 2pt and 4pt correlation functions
verboseMode = 1;
guessUpdateMode = 0;    % Very detailed: shows changes to any paramters
diagnoseMode = 1;       % Shows relative contribution of the 2pt TCF to the histogram at the end of each chi-squared calculation

clockMode = 0;          % Times various features of code.
saveModeDefault = 1;
plotMode = 1;           % Makes plots.
useFigurePosnMode = 1;

plot_Gen1_memberGuessesMode = 0; % Figure(1-2-3): Plots Histograms, C2, C4 of Gen1 (SLOW!)
plot_Gen1ChiSquaredMode = 1;     % Figure(4): Adds a point to the chi-sq vs member plot (plot4) (SLOW); Plot statistics from the first generation (informative)

plot_GenerationalProgressMode = 1; % Figure 4: Gray circles on Fig(4); Plot a generational update with best fit from the current generation
plot_paramater_histogram_mode = 0; % Figure 5

plot_model_GenerationalUpdatesMode = 1;       % Figure 6: Makes a model of the system with rates between states EACH GENERATIOON
plot_GenerationalChiSquaredHistogramMode = 0; % Figure 7: Histogram of chi-squared values
pauseBetweenGenerationMode = 0;               % Requires manual button press to continue

% At the end of the program
showBestFitMode = 1;  % Replaces current guess with best guess (end of code)

fitHistMode = 1;  fitHistData_mode = 1; % If 0 you will fit to the histogram fit
fitC2Mode = 1;
fitC4Mode = 1;

% Choose Weighting Amount
weightingFactor_FREThist = 1000;  % Weighting for FRET hist comparison
weightingFactor_C2 = 1;
weightingFactor_C4_t0 = 1000;

if sum([fitHistMode,fitC2Mode,fitC4Mode]) == 0
    error('No surfaces to optimize to! Pick at least one.')
end

%**************************************************************************
% REMOVE: ONLY HAVE NUMERICAL
% useAnalyticalAlgorithmsMode = 1;
% if useAnalyticalAlgorithmsMode == 1
%     if verboseMode == 1
%         disp('     ***using ANALYTICAL algorithms');
%     end
% else
%     if verboseMode == 1
%         disp('     ***using NUMERICAL algorithms');
%         disp('           ==> Loading the conditional probabilities as a function of rates');   % Change this to a function of eigenvector components?
%     end
%     load('filename')
% end

%**************************************************************************
ContributionArray = zeros(1,3);

if verboseMode == 1
    fprintf('       Each generation will have %d/%d of its members reproduce(%d%%)\r',Nreproduce,NmembersInitPop, percentReproduce);
    fprintf('       Each generation will see a total of %d mutations. (%d per individual on average).\r', Nmutations,Nmutations/NmembersInitPop);
    fprintf('       Each generation will have the top %d percent unmutated (%d/%d individuals).', percentToKeep,membersToKeep, NmembersInitPop);
end
%%
%//////////////////////////////////////////////////////////////////////////
% PART 1: Model Specific Paramaters
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% EIGHT STATE  MODEL HAS 3 LOOPS (DNA only, monomer, dimer)
% 1 <--> 2 <--> 3 <-->1                         (cyclic Dna-only)
%        2 <--> 5 <--> 6 <--> 3                 (monomer - productive path)
%        2 <--> 4             3 <-->7           (monomer - unproductive
%                      6 <--> 8 <--> 9 <--> 6   (dimer)
% The 3 state where the unbound DNA is fully extended has been omitted with
% because we predict that it will show up very rarely  in the system
% contaiing protein.

sigma_A1 = 0.15;
sigma_A2 = 0.1;
sigma_A3 = 0.1;
sigma_A4 = 0.1;
sigma_A5 = 0.1;
sigma_A6 = 0.1;
sigma_A7 = 0.1;
sigma_A8 = 0.1;
% sigma_A9 = 0.1;

sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4; sigma_A5; sigma_A6; sigma_A7; sigma_A8];

% Need to set bounds for each time constant
% Set DNA-only states from data fitting of 3state model
t12_bounds = [1e-6,3000e-6];        %Paramater #1 is high--> med
t13_bounds = [100e-6,500e-3];       %Paramater #2 is high --> low
t21_bounds = [1e-6,1e-3];           %Paramater #3 is med --> high
t23_bounds = [1e-6,10e-3];          %Paramater #4 is med --> low
t31_bounds = [10e-6,10e-3];         %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
% *t32 wll be determined by the others in loop

% Protein binding time constants:
t24_bounds = [1e-6, 1e-1];
t42_bounds = [1e-6, 1e-1];
t25_bounds = [1e-6, 1e-1];
t52_bounds = [1e-6, 1e-1];
t56_bounds = [1e-6, 1e-1];
% * t65_bounds will be determined by others in loop
t63_bounds = [1e-6, 1e-1];
t36_bounds = [1e-6, 1e-1];
t37_bounds = [1e-6, 1e-1];
t73_bounds = [1e-6, 1e-1];
t68_bounds = [1e-6, 1e-1];
t86_bounds = [1e-6, 1e-1];
% t89_bounds = [1e-6, 1e-1];
% t98_bounds = [1e-6, 1e-1];
% % * t96_bounds will be determined by others in loop
% t69_bounds = [1e-6, 1e-1];

% Define bounds of FRET parameters
A1_bounds = [0.65,1];   % Compact       (higest FRET)            %Paramater #6 % HIGH fret State
A2_bounds = [0.3,0.65]; % Not extended  (Intermediate FRET)      %Paramater #7 % Med FRET state
A3_bounds = [0,0.3];    % Extended      (lowest FRET)            %Paramater #8 %Low FRET state
A4_bounds = [0.3,0.65]; % Not extended
A5_bounds = [0.3,0.65]; % Not extended
A6_bounds = [0,0.3];    % Extended
A7_bounds = [0,0.3];    % Extended
A8_bounds = [0,0.3];    % Extended
% A9_bounds = [0,0.3]

boundsArray = [t12_bounds; t13_bounds; t21_bounds; t23_bounds; t31_bounds;...
    t24_bounds; t42_bounds; t25_bounds; t52_bounds; t56_bounds; t63_bounds;...
    t36_bounds; t37_bounds; t73_bounds; t68_bounds; t86_bounds;...
    A1_bounds; A2_bounds; A3_bounds; A4_bounds; A5_bounds;A6_bounds; A7_bounds;...
    A8_bounds;];

% boundsArray entries removed for 8 state (put back  for 9state):
% t89_bounds; t98_bounds; t69_bounds; A9_bounds;

Nparam = length(boundsArray);
maxMutationCountsArray = zeros(1, Nparam);
minMutationCountsArray = zeros(1, Nparam);

% Maximum mutation value
Max_mut_factor = 0.1;
t12_mutate = Max_mut_factor * t12_bounds(2)/2;%100e-6;
t13_mutate = Max_mut_factor * t13_bounds(2)/2;%100e-6;
t21_mutate = Max_mut_factor * t21_bounds(2)/2;%500e-6;
t23_mutate = Max_mut_factor * t23_bounds(2)/2;%10e-6;
t31_mutate = Max_mut_factor * t31_bounds(2)/2;%100e-6;
% t32 -> detailed balance
t24_mutate = Max_mut_factor * t24_bounds(2)/2;
t42_mutate = Max_mut_factor * t42_bounds(2)/2;
t25_mutate = Max_mut_factor * t25_bounds(2)/2;
t52_mutate = Max_mut_factor * t52_bounds(2)/2;
t56_mutate = Max_mut_factor * t56_bounds(2)/2; 
% t65  -> detailed balance
t63_mutate = Max_mut_factor * t63_bounds(2)/2;
t36_mutate = Max_mut_factor * t36_bounds(2)/2;
t37_mutate = Max_mut_factor * t37_bounds(2)/2; 
t73_mutate = Max_mut_factor * t73_bounds(2)/2;
t68_mutate = Max_mut_factor * t68_bounds(2)/2;
t86_mutate = Max_mut_factor * t86_bounds(2)/2; 
% t89_mutate = Max_mut_factor * t89_bounds(2)/2;
% t98_mutate = Max_mut_factor * t98_bounds(2)/2;
% % * t96 -> detailed balance
% t69_mutate = Max_mut_factor * t69_bounds(2)/2; 

A1_mutate = (A1_bounds(2) - A1_bounds(1))*0.1;%0.05;
A2_mutate = (A2_bounds(2) - A2_bounds(1))*0.1;%0.05;0.05;
A3_mutate = (A3_bounds(2) - A3_bounds(1))*0.1;%0.05;0.05;
A4_mutate = (A4_bounds(2) - A4_bounds(1))*0.1;
A5_mutate = (A5_bounds(2) - A5_bounds(1))*0.1;
A6_mutate = (A6_bounds(2) - A6_bounds(1))*0.1;
A7_mutate = (A7_bounds(2) - A7_bounds(1))*0.1;
A8_mutate = (A8_bounds(2) - A8_bounds(1))*0.1;
% A9_mutate = (A9_bounds(2) - A9_bounds(1))*0.1;

maxMutationArray = [t12_mutate;t13_mutate;t21_mutate;t23_mutate;t31_mutate;...
    t24_mutate; t42_mutate; t25_mutate; t52_mutate; t56_mutate; ...
    t63_mutate; t36_mutate; t37_mutate; t73_mutate; ...
    t68_mutate; t86_mutate; ...
    A1_mutate;A2_mutate;A3_mutate; A4_mutate; A5_mutate; A6_mutate; ...
    A7_mutate; A8_mutate;];% A9_mutate];  t89_mutate; t98_mutate; t69_mutate;

%Start in the single molecule folder (smData_Processed): comp specific
[computer_terminal_str, terminalID] = computerMode(pwd);

NconstructFolderNames = numel(constructFolderNames);

for construct_idx = 1:NconstructFolderNames
    constructName = char(constructFolderNames(construct_idx));
    disp(['     Construct' num2str(construct_idx) '/' num2str(NconstructFolderNames) ': ' constructName]);
    
    if exist('constructName', 'var') ~= 1
        error('Not sure what construct to analyze');
    end
    
    % Find GenAlg files and load target data:
    [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix] = fileLocator(terminalID, constructName, modelName, protein_str);

    % Load file path for data plots
    [histogram_FilePath, C2_FilePath, C4_FilePath] = loadFilePath_dataPlots(genAlgDataLocationPrefix, protein_str);   
  
  %--------------------------------------------------------------------------
    % Set up the figure positions once and done
    %--------------------------------------------------------------------------
    if plotMode == 1
        clear figPosn
        if useFigurePosnMode == 1
            figPosn = setFigPosn(terminalID);
            if fitHistMode == 1
                figure(1);
                clf;
                set(gcf,'Name','FRET Histogram');
                set(gcf,'Color','w');
                set(gcf,'Position',figPosn(:,1));
            end
            %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if fitC2Mode == 1
                figure(2);
                clf;
                set(gcf,'Name','Two-Point TCF: C2');
                set(gcf,'Position',figPosn(:,2));
            end
            %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if fitC4Mode == 1
                figure(3)
                clf;
                set(gcf,'Name','Four-Point TCF: C4');
                set(gcf,'Position',figPosn(:,3));
            end
            %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if plot_model_GenerationalUpdatesMode == 1
                figure(4);
                clf;
                set(gcf,'Name','Fitting Progress');
                set(gcf,'Position',figPosn(:,4));
            end
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        if plot_paramater_histogram_mode == 1
            figure(5);
            clf;
            set(gcf,'Name',' Distribution of initial paramaters : Gen 1');
            set(gcf,'Position',figPosn(:,5));
            set(gcf,'Color','w');
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        % Fig 6: Model of the network with rates and FRET values
        if plot_model_GenerationalUpdatesMode == 1
            figure(6);
            clf;
            set(gcf,'Color','w');
            set(gcf,'Name','Model: 3state123 cyclical');
            if useFigurePosnMode == 1
                set(gcf,'Position',figPosn(:,6));
            end
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        if plot_GenerationalChiSquaredHistogramMode == 1
            figure(7);
            clf;
            set(gcf,'Name','Chi-squared values guesses');
            set(gcf,'Position',figPosn(:,7));
            set(gcf,'Color','w');
        end
        
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    end    
%% Part 2: Load data     
    %--------------------------------------------------------------------------
    % (1) Optimization Target #1: 1-D FRET HISTOGRAM (PART 2: Load the target data)
    %--------------------------------------------------------------------------
    % if fitHistMode == 1
    %     load(histogram_FileName,'xData','yData','xFit','yFit');
    load(histogram_FilePath,'xData','yData');
    if fitHistData_mode == 1
        FRET_bins = xData;
        targetHistogram = yData;
    elseif fitHistData_mode == 0
        FRET_bins = xFit;
        targetHistogram = yFit;
    end
    targetHistogram = targetHistogram./sum(targetHistogram);
    
    if plotMode == 1
        figure(1);
        set(gcf,'Name','FRET Histogram');
        set(gcf,'Color','w');
        
        data_hist_Plot = plot(FRET_bins,targetHistogram);
        data_hist_Plot.LineWidth = 2;
        data_hist_Plot.Color = 'blue';
%         x.DisplayName = 'Data';
        
        xlabel('FRET Efficiency','FontSize',14);
        ylabel('Frequency','FontSize',14);
        
        [sample_description, ~] = sample_descriptionGetter();
        title_str = ['Experimental vs Simulated Histograms' ...
            10 sample_description];
        title(title_str,'fontsize',14);
        
        
        
        lgd = legend('show');
        lgd.Location = 'northwest';
        lgd.FontSize = 12;
        drawnow();
        hold on;
    end
    
    % end
    
    %--------------------------------------------------------------------------
    % (2) Optimization Target #2: 2-point TCF (C2)           (20 uSec and on)  (PART 2: Load the target data)
    %--------------------------------------------------------------------------
    time  =  time_sim;
    if fitC2Mode == 1
%         load(C2_FilePath,'time','yData','y','yoff');
        load(C2_FilePath,'tau_Sec_array','C2_array','yoff');
        fitC2Data_mode = 1; % If 0 you will fit to the histogram fit
        if fitC2Data_mode == 1
            C2_exp_x = tau_Sec_array;
            C2_time = C2_exp_x;
            C2_exp_y = C2_array;
        elseif fitC2Data_mode == 0
            C2_exp_x = tau_Sec_array;
            C2_time = C2_exp_x;
%             C2_exp_y = y;
        end
        if normalizeMode == 1
            C2_exp_y = C2_exp_y./C2_exp_y(1);
        end
        weightC2func = 1./(sqrt(C2_exp_x));% 2-pt. TCF weighting function
        if plotMode == 1
            figure(2)
            set(gcf,'Color','w');
            hold on;
            plot(C2_exp_x,C2_exp_y,'b.','MarkerSize',10,'DisplayName','C^{(2)}(\tau) Data');
            
            [sample_description, ~] = sample_descriptionGetter();
            title_str = ['Experimental vs Simulated C2' ...
                10 sample_description];
            title(title_str,'fontsize',14);
            
            xlabel('Time (sec)','FontSize',14);
            ylabel('C^{(2)}(\tau)','FontSize',14);
            set(gca,'xscale','log');
            drawnow();
        end
    end
    
    %---------------------------------------------------------------------------------------------
    % (3) Optimization Target #3: 4-point TCF (C4)    (10 uSec and on) (PART 2: Load the target data)
    %---------------------------------------------------------------------------------------------
    if fitC4Mode == 1
        %         load(C4_FilePath,'FourPtTCF_avg','tau1arrayUsec','tau3arrayUsec','tau2ValUsec');
        %  C4_tau2eq0_exp = FourPtTCF_avg;
        %         C4_tau1range = tau1arrayUsec*1e-6;
        %         C4_tau3range = tau3arrayUsec*1e-6';
        
        load(C4_FilePath,'fourptTCF','tau1arraySec','tau3arraySec','tau2ValSec');
        C4_tau2eq0_exp = fourptTCF;  % was  called C4 before
        C4_tau1range = tau1arraySec;
        C4_tau3range = tau3arraySec';
        C4_time = tau1arraySec;
        tau2 = tau2ValSec;
        
        load(C2_FilePath,'yoff'); 
        zoff = yoff*yoff;
        if normalizeMode == 1
            C4_tau2eq0_exp = C4_tau2eq0_exp./C4_tau2eq0_exp(1,1);
        end
        
        % 4-pt. TCF weighting function
        wC4func = 1./sqrt(C4_tau1range).*(1./(sqrt(C4_tau1range)));
        
        if plotMode == 1
            figure(3)
            clf;
            set(gcf,'Name','Four-Point TCF: C^{(4)}');
            set(gcf,'Color','w');
            
            C4dataPlot = mesh(C4_tau1range,C4_tau3range,C4_tau2eq0_exp);
            C4dataPlot.DisplayName = 'C^{(4)} Data';
            
            xlabel('\tau_1 (sec)','FontSize',14);
            ylabel('\tau_3 (sec)','FontSize',14);
            ylabel('C^{(4)}','FontSize',14);
            %         title('Experimental vs Simulated C4','FontSize',14);
            
            [sample_description, save_prefix] = sample_descriptionGetter();
            title_str = ['C^{(4)}(\tau_1, \tau_2 = ' num2str(tau2ValSec*1e6) '\musec, \tau_3)' ...
                10 sample_description];
            title(title_str,'fontsize',14);
            xlabel('\tau_1 (sec)','fontsize',14);
            ylabel('\tau_3 (sec)','fontsize',14);
            
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
    
    if verboseMode == 1
        disp('     Part1: Done Loading and (optionally) plotting the data.');
    end
    % disp('Press enter to continue');
    % pause();
    
    folderInName = [programName '_output' filesep() 'lowestChiSquare'];
    finName =  'BestFitResults.mat';
    finPath = [folderInName filesep() finName];
    if exist([programName '_output' filesep() 'lowestChiSquare'],'dir') == 7
        load(finPath,'iter');
        if forceMoreIterationsMode == 1
            if iter > maxIterations
                maxIterations = iter + NumberToExtend;
            end
        end
    else
        iter = 0;
    end
    
    while iter < maxIterations
        saveMode = saveModeDefault; %Saves the output of the program to a folder
        
        iter = iter + 1;
        disp(['iter = ' num2str(iter) ]);
        
%% Part 3: Gen1
        %//////////////////////////////////////////////////////////////////////////
        % PART 3: Make the initial generation of guesses           (Part3: Gen1)
        %//////////////////////////////////////////////////////////////////////////
        %Give each member of the population some value between the LB and UB
        population = rand(NmembersInitPop,Nparam);
        for param_idx = 1:Nparam
            %To pick a random number in the interval of LB to UB:
            % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
            population(:,param_idx) = boundsArray(param_idx,1) + population(:,param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
        end
        
        newpopulation = zeros(size(population));
        pop_chisquared_array = zeros(NmembersInitPop,1);
        guess = zeros(maxGenerations,9);
        genNum = 1;
        tic
%         parfor pop_idx = 1:NmembersInitPop
        for pop_idx = 1:NmembersInitPop
            disp('member')
            disp(pop_idx)
            t12 = population(pop_idx,1);
            t13 = population(pop_idx,2);
            t21 = population(pop_idx,3);
            t23 = population(pop_idx,4);
            t31 = population(pop_idx,5);
            t24 = population(pop_idx,6);
            t42 = population(pop_idx,7);
            t25 = population(pop_idx,8);
            t52 = population(pop_idx,9);
            t56 = population(pop_idx,10);
            t65 = population(pop_idx,11);
            % * t63_bounds will be determined by others in loop
            t36 = population(pop_idx,12);
            t37 = population(pop_idx,13);
            t73 = population(pop_idx,14);
            t68 = population(pop_idx,15);
            t86 = population(pop_idx,16);
%             t89 = population(pop_idx,17);
%             t98 = population(pop_idx,18);
            % * t96_bounds will be determined by others in loop
%             t69 = population(pop_idx,20);

%             % Loop conditions
%             t32 = (t12 * t23 * t31)/(t13 * t21);
%             t63 = (t23 * t36 * t65* t52)/(t32 * t25 * t56);
%            
                                            % Need these indexes for 9 state
            A1 = population(pop_idx,17);        %21);
            A2 = population(pop_idx,18);        %22);
            A3 = population(pop_idx,19);        %23);
            A4 = population(pop_idx,20);        %24);
            A5 = population(pop_idx,21);        %25);
            A6 = population(pop_idx,22);        %26);
            A7 = population(pop_idx,23);        %27);
            A8 = population(pop_idx,24);        %28);
%             A9 = population(pop_idx,29);

          % Define rates, kij, from the tij's
          k12 = 1/t12;
          k13 = 1/t13;
          k21 = 1/t21;
          k23 = 1/t23;
          k31 = 1/t31;
          k24 = 1/t24;
          k42 = 1/t42;
          k25 = 1/t25;
          k52 = 1/t52;
          k56 = 1/t56;
          k65 = 1/t65;
          k36 = 1/t36;
          k37 = 1/t37;
          k73 = 1/t73;
          k68 = 1/t68;
          k86 = 1/t86;
          %           k89 = 1/t89;
          %           k98 = 1/t98;
          %           t96 = 1/t96;
          %           t69 = 1/t69;
          
          % Loop conditions
          k32 = (k12 * k23 * k31)/(k13 * k21);
          k63 = (k23 * k36 * k65* k52)/(k32 * k25 * k56);
          t32 = 1/k32;
          t63 = 1/k63;

         A = [A1;A2;A3;A4;A5;A6;A7;A8];
         rates = [k12,k13,k21,k23,k32,k31,k24,k42,k25,k52,k56,k65,k63,k36,k37,k73,k68,k86];

         % Define K matrix
         % 8 state model:
         K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0;...
             k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0;...
             k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0;...
             0, k24, 0, -k42, 0, 0, 0, 0;...
             0, k25, 0, 0, -(k52+k56), k65, 0, 0;...
             0, 0, k36, 0, k56, -(k65+k63+k68), 0, k86;...
             0, 0, k37, 0, 0, 0, -k73, 0;...
             0, 0, 0, 0, 0, k68, 0, -k86;];
         
%          % 9 state model:
%          K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0, 0;...
%              k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0, 0;...
%              k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0, 0;...
%              0, k24, 0, -k42, 0, 0, 0, 0, 0;...
%              0, k25, 0, 0, -(k52+k56), k65, 0, 0, 0;...
%              0, 0, k36, 0, k56, -(k65+k63+k68+k69), 0, k86, k96;...
%              0, 0, k37, 0, 0, 0, -k73, 0, 0;...
%              0, 0, 0, 0, 0, k68, 0, -(k86+k89),k98;...
%              0, 0, 0, 0, 0, k69, 0, k89, -(k96+k98)];
%          
         
         % Define the conditional probability matrix:
         %          syms  t
         %          p = k2P(K,t);
         time =  time_sim;
         [P, ~, p, time] = k2P(K,time);
         
         N = length(K);
 %%        
          if plotMode == 1
                %--------------------------------------------------------------------------
                % gen1: Plot the array of initial guesses             (PART 3: Gen1)
                %--------------------------------------------------------------------------
                if plot_Gen1_memberGuessesMode == 1
                    %--------------------------------------------------------------------------
                    % (1) Optimization Target #1: 1-D FRET HISTOGRAM  (PART 3: Gen1)
                    %--------------------------------------------------------------------------
                    if fitHistMode == 1
                        figure(1);
                        [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate(P, A, sigma_A);
                        
                        if plotMode == 1
                            figure(1);
                            clf;
                            
                            set(gcf,'Color','w');
                            set(gcf,'Name','FRET Histogram');
                            lineStyle = char('r--','g--','b--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                            
                            FRET_bins = linspace(0,1,100);
                            
                            if exist('data_hist_Plot','var') == 1
                                delete(data_hist_Plot)
%                                 disp('data_hist_Plot deleted')
                            end
                            % Replot data Histogram 
                            data_hist_Plot = plot(FRET_bins,targetHistogram);
                            data_hist_Plot.LineWidth = 2;
                            data_hist_Plot.Color = 'blue';
                            x.DisplayName = 'Data';
                            
                            xlabel('FRET Efficiency','FontSize',14);
                            ylabel('Frequency','FontSize',14);
                            
                            [sample_description, ~] = sample_descriptionGetter();
                            title_str = ['Experimental vs Simulated Histograms' ...
                                10 sample_description];
                            title(title_str,'fontsize',14);
                            
                            hold on;

                            % Clear plots from previous run
                            if exist('histPlot','var') == 1
                                delete(histPlot)
                                disp('histPlot deleted')
                            end
                            
                            % Plot histogram of each state
                            for i = 1:N
                                histPlot = plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
                                    lineStyle(i,:),'LineWidth',1,'DisplayName',Peq_names(i,:));
                                lgd = legend('show');
                                % lgd.Location = 'northwest';
                                lgd.FontSize = 14;
                                hold on
                            end
                            
                            set(gca, 'FontSize', 14);
                            

                            if exist('histPlotTot','var') == 1
                                delete(histPlotTot)
                            end
                            
                            hold on
                            
                            % Plot overall histogram
                            histPlotTot = plot(FRET_bins, hist_sim,...
                                lineStyle(i),'LineWidth',1,'DisplayName','Total Fit');
                            lgd_tot = legend('show');
                            
                        end
                    end
                    %---------------------------------------------------------------------------------
                    % (2) Optimization Target #2: 2-point TCF (C2)  (20 uSec and on)  (PART 3: GEN1)
                    %---------------------------------------------------------------------------------
                    if fitC2Mode == 1
                        [time, C2, ~] = P2C(P, K, time);
                        if plotMode == 1
                            figure(2)
                            set(gcf,'Color','w');
                            
                            if exist('C2_plot','var')  == 1
                                delete(C2_plot)
                            end
                            
                            C2_plot = plot(time,C2);
                            title_str = ['Two point time correlation function'];
                            title(title_str,'FontSize',18);
                            xlabel('\tau (sec)','fontsize',16);
                            ylabel('C^{(2)}(\tau)','fontsize',16);
                            set(gca,'yscale','linear');
                            set(gca,'xscale','log');
                            set(gca,'FontSize',14);
                            grid on
                            axis tight;
                            
                            drawnow();
                        end
                        
                    end
                    %-----------------------------------------------------------------------------
                    % (3) Optimization Target #3: 4-point TCF (C4) (20 uSec and on) (PART 3: GEN1)
                    %------------------------------------------------------------------------------
                    if fitC4Mode == 1
%                         if exist('C4_plot','var')  == 1
%                             delete(C2_plot)
%                         end
                        [time,~,C4] = P2C(P, K, time);
                        
                        if plotMode == 1
                            figure(3)
                            
                            if exist('C4_plot','var')  == 1
                                delete(C4_plot)
                            end
                            
                            set(gcf,'Color','w');
                            hold on;
                            C4_plot = surf(time,time,C4);
                            title_str = ['Four point time correlation function'];
                            title(title_str,'FontSize',18);
                            xlabel('\tau_1 (sec)','fontsize',16);
                            ylabel('\tau_3 (sec)','fontsize',16);
                            zlabel('C^{(4)}(\tau)','fontsize',16);
                            set(gca,'yscale','log');
                            set(gca,'xscale','log');
                            set(gca,'FontSize',14);
                            grid on
                            axis tight;
                            
                            drawnow();
                        end
                    end
                end % End of plotting guesses for gen1
          end % end of plot mode
          %--------------------------------------------------------------------------
          % Asses the chisquared of the guess                 (PART 3: Gen1)
          %--------------------------------------------------------------------------
          [chisquared,chisquared_array] = chiSqCalc(rates, A, P, sigma_A, C2_time, C4_time);
          
          pop_chisquared_array(pop_idx) = chisquared;
          
        end%End of "pop_idx: Loop
        if verboseMode == 1
            disp('     Part2 : Done simulating the first generation. Selecting the best guess to plot.');
        end
        elapsedTime = toc;
        if clockMode == 1
            disp(['     It took ' num2str(elapsedTime) ' seconds to compute the '...
                'RMS of ' num2str(NmembersInitPop) ' guesses.(' num2str(elapsedTime/NmembersInitPop) ' sec/guess)']);
        end
        
        %------------------------------------------------------------------------------
        % (4) Plot the chisquared value of each member of the population (PART 3: Gen1)
        %-------------------------------------------------------------------------------
        if plotMode == 1
            if plot_Gen1ChiSquaredMode == 1
                figure(4);
                  
                set(gcf,'Color','w');
                set(gca,'yscale','log');
                for pop_idx = 1:length(pop_chisquared_array)
                    semilogy(pop_idx,pop_chisquared_array(pop_idx),'--gs',...
                        'LineWidth',2,...
                        'MarkerSize',10,...
                        'MarkerEdgeColor','b',...
                        'MarkerFaceColor',[0.5,0.5,0.5]);
                    hold on;
                end
                hold on;
                axis tight;
                xlim([1,inf]);
                
                xlabel('Member # (Gen1)','FontSize',14);
                ylabel('RMS Fitness (\chi^2)','FontSize',14);
                title('\chi^2 of Gen 1','FontSize',14)
                drawnow();
            end
        end
%%                
        %-------------------------------------------------------------------------------
        % (7) Make a histogram of the chisquared values of each member   (PART 3: Gen1)
        %--------------------------------------------------------------------------------
        if plotMode == 1
            if plot_GenerationalChiSquaredHistogramMode == 1
                figure(7)
                set(gcf,'Name','Chi-squared values of gen1 guesses');
                
                initialHist = histogram(pop_chisquared_array);
                xlabel('\chi^2 Value');
                ylabel('Frequency');
                title('Distribution of Gen1 guesses chi-squared');
            end
        end
        
        %--------------------------------------------------------------------------
        % Sort Generation 1 as a function of chisquare (PART 3: Gen1)
        %--------------------------------------------------------------------------
        % Sort from Lowest rms to the highest rms (ascending order)
        
        [pop_chisquare_sorted,index_arr] = sort(pop_chisquared_array);%First in array have the lowest chisquare value
        
        population = population(index_arr,:);
        
        guess(genNum,1:Nparam) = population(1,:); % Current best guess ****************************************************
        
        Best_chisquared = pop_chisquare_sorted(1);
        guess(genNum,Nparam+1) = Best_chisquared; % Current best guess' fit value
        
        genNum_array = genNum;
        chisquared_array = Best_chisquared;
                
        % Display best guess to screen
        
        t12 = guess(genNum,1);
        t13 = guess(genNum,2);
        t21 = guess(genNum,3);
        t23 = guess(genNum,4);
        t31 = guess(genNum,5);
        % * t32 will be determined by others in loop
        t24 = guess(genNum,6);
        t42 = guess(genNum,7);
        t25 = guess(genNum,8);
        t52 = guess(genNum,9);
        t56 = guess(genNum,10);
        t65 = guess(genNum,11);
        % * t63 will be determined by others in loop
        t36 = guess(genNum,12);
        t37 = guess(genNum,13);
        t73 = guess(genNum,14);
        t68 = guess(genNum,15);
        t86 = guess(genNum,16);
        
        % Loop condition
%         t32 = (t12 * t23 * t31)/(t13 * t21);
%         t63 = (t23 * t36 * t65* t52)/(t32 * t25 * t56);
        
        A1 = guess(genNum,17);        %21);
        A2 = guess(genNum,18);        %22);
        A3 = guess(genNum,19);        %23);
        A4 = guess(genNum,20);        %24);
        A5 = guess(genNum,21);        %25);
        A6 = guess(genNum,22);        %26);
        A7 = guess(genNum,23);        %27);
        A8 = guess(genNum,24);
        
        chisquared = guess(genNum,25);
        
        % Redefine rates, kij, from the tij's
          k12 = 1/t12;
          k13 = 1/t13;
          k21 = 1/t21;
          k23 = 1/t23;
          k32 = 1/t32;
          k31 = 1/t31;
          k24 = 1/t24;
          k42 = 1/t42;
          k25 = 1/t25;
          k52 = 1/t52;
          k56 = 1/t56;
          k65 = 1/t65;
          k63 = 1/t63;
          k36 = 1/t36;
          k37 = 1/t37;
          k73 = 1/t73;
          k68 = 1/t68;
          k86 = 1/t86;
%           k89 = 1/t89;
%           k98 = 1/t98;
%           t96 = 1/t96;
%           t69 = 1/t69;

% Recalculate loop conditions
          k32 = (k12 * k23 * k31)/(k13 * k21);
          t32 = 1/k32;
          k63 = (k23 * k36 * k65* k52)/(k32 * k25 * k56);
          t63 = 1/k63;
          
        A = [A1;A2;A3;A4;A5;A6;A7;A8];
        rates = [k12,k13,k21,k23,k32,k31,k24,k42,k25,k52,k56,k65,k63,k36,k37,k73,k68,k86];

        if verboseMode == 1 
            fprintf(['Best fit from initial generation was member #%d:\n '...
                't12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
                't24 = %f, t42 = %f, t25 = %f, t52 = %f, t56 = %f, t65 = %f'...
                't63 = %f, t36 = %f,t37 = %f, t73 = %f,t68 = %f, t86 = %f'...
                '\n A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f,'...
                'A6 = %f, A7 = %f, A8 = %f,Best_chisquared = %f\r\n'],...
                index_arr(1),t12,t13,t21,t23,t31,t32,...
                t24, t42, t25, t52, t56, t65, t63, t36,...
                t37, t73, t68, t86,A1, A2, A3,A4, A5, A6, A7, A8,...
                Best_chisquared);
        end
        % Define K matrix
         % 8 state model:
         K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0;...
             k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0;...
             k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0;...
             0, k24, 0, -k42, 0, 0, 0, 0;...
             0, k25, 0, 0, -(k52+k56), k65, 0, 0;...
             0, 0, k36, 0, k56, -(k65+k63+k68), 0, k86;...
             0, 0, k37, 0, 0, 0, -k73, 0;...
             0, 0, 0, 0, 0, k68, 0, -k86;];
         
%          % 9 state model:
%          K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0, 0;...
%              k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0, 0;...
%              k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0, 0;...
%              0, k24, 0, -k42, 0, 0, 0, 0, 0;...
%              0, k25, 0, 0, -(k52+k56), k65, 0, 0, 0;...
%              0, 0, k36, 0, k56, -(k65+k63+k68+k69), 0, k86, k96;...
%              0, 0, k37, 0, 0, 0, -k73, 0, 0;...
%              0, 0, 0, 0, 0, k68, 0, -(k86+k89),k98;...
%              0, 0, 0, 0, 0, k69, 0, k89, -(k96+k98)];

            % Deine our new  P using this best guess K:
            [P, ~, ~, time] = k2P(K,time);
%%
        %--------------------------------------------------------------------------
        % (1) Make a histogram which corresponds to the best guess (PART 3: Gen1)
        %--------------------------------------------------------------------------
        if plotMode == 1
            %--------------------------------------------------------------------------
            % (1) Optimization Target #1: 1-D FRET HISTOGRAM  (PART 3: Gen1)
            %--------------------------------------------------------------------------
            if fitHistMode == 1
                figure(1)
                hold on;
                [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate(P, A, sigma_A);
                
                if plotMode == 1
                            figure(1);
                            clf;
                            
                            set(gcf,'Color','w');
                            set(gcf,'Name','FRET Histogram');
                            lineStyle = char('r--','g--','b--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                            
                            FRET_bins = linspace(0,1,100);
                            
                            if exist('data_hist_Plot','var') == 1
                                delete(data_hist_Plot)
%                                 disp('data_hist_Plot deleted')
                            end
                            % Replot data Histogram 
                            data_hist_Plot = plot(FRET_bins,targetHistogram);
                            data_hist_Plot.LineWidth = 2;
                            data_hist_Plot.Color = 'blue';
%                             x.DisplayName = 'Data';
                            
                            xlabel('FRET Efficiency','FontSize',14);
                            ylabel('Frequency','FontSize',14);
                            
                            [sample_description, ~] = sample_descriptionGetter();
                            title_str = ['Experimental vs Simulated Histograms' ...
                                10 sample_description];
                            title(title_str,'fontsize',14);
                            
                            hold on;

                            % Clear plots from previous run
                            if exist('histPlot','var') == 1
                                delete(histPlot)
                                disp('histPlot deleted')
                            end
                            
                            % Plot histogram of each state
                            for i = 1:N
                                histPlot = plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
                                    lineStyle(i,:),'LineWidth',1,'DisplayName',Peq_names(i,:));
                                lgd = legend('show');
                                % lgd.Location = 'northwest';
                                lgd.FontSize = 14;
                                hold on
                            end
                            set(gca, 'FontSize', 14);
                            if exist('histPlotTot','var') == 1
                                delete(histPlotTot)
                            end
                            hold on
                            % Plot overall histogram
                            histPlotTot = plot(FRET_bins, hist_sim,...
                                lineStyle(i),'LineWidth',1,'DisplayName','Total Fit');
                            lgd_tot = legend('show');
                            
                end
            end
                %---------------------------------------------------------------------------------
                % (2) Optimization Target #2: 2-point TCF (C2)  (20 uSec and on)  (PART 3: GEN1)
                %---------------------------------------------------------------------------------
                if fitC2Mode == 1
                    [time, C2, ~] = P2C(P, K, time);
                    if plotMode == 1
                        figure(2)
                        set(gcf,'Color','w');
                        
                        if exist('C2_plot','var')  == 1
                            delete(C2_plot)
                        end
                        
                        C2_plot = plot(time,C2);
                        title_str = ['Two point time correlation function'];
                        title(title_str,'FontSize',18);
                        xlabel('\tau (sec)','fontsize',16);
                        ylabel('C^{(2)}(\tau)','fontsize',16);
                        set(gca,'yscale','linear');
                        set(gca,'xscale','log');
                        set(gca,'FontSize',14);
                        grid on
                        axis tight;
                        
                        drawnow();
                    end  
                end
                %-----------------------------------------------------------------------------
                % (3) Optimization Target #3: 4-point TCF (C4) (20 uSec and on) (PART 3: GEN1)
                %------------------------------------------------------------------------------
                if fitC4Mode == 1
                    %                         if exist('C4_plot','var')  == 1
                    %                             delete(C2_plot)
                    %                         end
                    [time,~,C4] = P2C(P, K, time);
                    
                    if plotMode == 1
                        figure(3)
                        
                        if exist('C4_plot','var')  == 1
                            delete(C4_plot)
                        end
                        
                        set(gcf,'Color','w');
                        hold on;
                        C4_plot = surf(time,time,C4);
                        title_str = ['Four point time correlation function'];
                        title(title_str,'FontSize',18);
                        xlabel('\tau_1 (sec)','fontsize',16);
                        ylabel('\tau_3 (sec)','fontsize',16);
                        zlabel('C^{(4)}(\tau)','fontsize',16);
                        set(gca,'yscale','log');
                        set(gca,'xscale','log');
                        set(gca,'FontSize',14);
                        grid on
                        axis tight;
                        
                        drawnow();
                    end
                end
                    %-------------------------------------------------------------------------------------------
                    % (4) Make a red circle around the best fit and begin optimization from there (PART 3: Gen1)
                    %-------------------------------------------------------------------------------------------
                    if plot_Gen1ChiSquaredMode == 1
                        figure(4);
                        plot(index_arr(1),chisquared,'--go','LineWidth',3,...
                            'MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','none')
                        y = [chisquared chisquared];
                        x = [index_arr(1) 1];
                        line1 = line(x,y);
                        line1.LineWidth = 2;
                        line1.Color = 'r';
                        line1.LineStyle = '--';
                    end
                    
                    %------------------------------------------------------------------------------
                    % (4) Keep track of the chisquared  for each member          (PART 3: Gen1)
                    %------------------------------------------------------------------------------
                    if plot_GenerationalProgressMode == 1
                        figure(4);
                        clf;
                        plot(genNum,chisquared,'bs',...
                            'MarkerSize',10,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor',[0.5,0.5,0.5]);
                        axis tight;
                        xlim([1,inf]);
                        set(gcf,'Color','w');
                        
                        xlabel('Generation Number','FontSize',14);
                        ylabel('RMS Fitness (\chi^2)','FontSize',14);
                        title('Fitness vs Generation Number','FontSize',14)
                        drawnow();
                    end
                    %--------------------------------------------------------------------------
                    % (5) Plot histograms of initial guesses for each paramater (PART 3: Gen1)
                    %--------------------------------------------------------------------------
                    if plot_paramater_histogram_mode == 1
                        figure(5)
                        param_strings = {'t12','t13','t21','t23','t31','t24','t42','t25','t52',...
                            't56','t65','t36','t37','t73','t68','t86','t32','t63',...
                            'A1','A2','A3','A4','A5','A6','A7','A8'};
                        maxCountArray = zeros(1,Nparam);
                        N_tij = Nparam - N; % N = number of states, N_tij is number of times *********************************
                        for param_idx = 1:Nparam
                            subplot(2,4,param_idx)
                            
                            if ismember(param_idx,1:N_tij)                      %Rate histograms
                                nbins = 50;
                                histogram(population(:,param_idx),nbins);
                                [N,edges] = histcounts(population(:,param_idx),nbins);
                                maxCountArray(param_idx) = max(N);
                                xlabel('Time (Sec)');
                            elseif ismember(param_idx,(N_tij+1):Nparam) == 1    %FRET Historgrams
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
                    end % if plot_parameter_histogram_mode
        end % if plotMode
        
        % error('Quitting before generation 2');
        if pauseBetweenGenerationMode == 1
            disp(['     Gen' num2str(genNum) ' Chi^2 = ' num2str(chisquared) ': Press Enter to proceed to the next generation']);
            pause();
        end
 %%
        %///////////////////////////////////////////////////////////////////////////////////////
        % PART 4: Refine the guesses by mutating certain paramaters randomly (PART 4: GenX)
        %///////////////////////////////////////////////////////////////////////////////////////
        %Bookmarking for keeping up with the number of trials
        Ntrials = 0;
        for genNum = 2:maxGenerations
            if clockMode == 1, tic; end
            % Mix genes of top rated individuals (population is sorted by chisquared: The lowest chisquared are first)
            %             newpop = [];
            %             newpopulation_temp = 0;
            for n = 1:NmembersInitPop - membersToKeep
                newpop = [];
                newpopulation_temp = 0;
                for i = 1:Nparam
                    newpopulation_temp = [population(ceil(Nreproduce*rand(1,1)),i)];
                    newpop = [newpop, newpopulation_temp];
                end
                newpopulation(n,:) = [newpop];
            end
%             for n = 1:NmembersInitPop - membersToKeep
%                 newpopulation(n,:) = [population(ceil(Nreproduce*rand(1,1)),1),...
%                     population(ceil(Nreproduce*rand(1,1)),2),...
%                     population(ceil(Nreproduce*rand(1,1)),3),...
%                     population(ceil(Nreproduce*rand(1,1)),4),...
%                     population(ceil(Nreproduce*rand(1,1)),5),...
%                     population(ceil(Nreproduce*rand(1,1)),6),...
%                     population(ceil(Nreproduce*rand(1,1)),7),...
%                     population(ceil(Nreproduce*rand(1,1)),8)];
%             end
            
            pop_idx = 0;
            while pop_idx <= membersToKeep
                %         Let the highest members be the best guesses from the first generation
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
        
        
        
        %--------------------------------------------------------------------------------
        %  After Mutations, Recalculate the chisquared values of each member (PART 4: GenX)
        %--------------------------------------------------------------------------------
        
        for pop_idx = 1:NmembersInitPop
            t12 = population(pop_idx,1);
            t13 = population(pop_idx,2);
            t21 = population(pop_idx,3);
            t23 = population(pop_idx,4);
            t31 = population(pop_idx,5);
            t24 = population(pop_idx,6);
            t42 = population(pop_idx,7);
            t25 = population(pop_idx,8);
            t52 = population(pop_idx,9);
            t56 = population(pop_idx,10);
            t65 = population(pop_idx,11);
            % * t63_bounds will be determined by others in loop
            t36 = population(pop_idx,12);
            t37 = population(pop_idx,13);
            t73 = population(pop_idx,14);
            t68 = population(pop_idx,15);
            t86 = population(pop_idx,16);
            %             t89 = population(pop_idx,17);
            %             t98 = population(pop_idx,18);
            % * t96_bounds will be determined by others in loop
            %             t69 = population(pop_idx,20);
            
            %             % Loop conditions
            %             t32 = (t12 * t23 * t31)/(t13 * t21);
            %             t63 = (t23 * t36 * t65* t52)/(t32 * t25 * t56);
            %
            % Need these indexes for 9 state
            A1 = population(pop_idx,17);        %21);
            A2 = population(pop_idx,18);        %22);
            A3 = population(pop_idx,19);        %23);
            A4 = population(pop_idx,20);        %24);
            A5 = population(pop_idx,21);        %25);
            A6 = population(pop_idx,22);        %26);
            A7 = population(pop_idx,23);        %27);
            A8 = population(pop_idx,24);        %28);
            %             A9 = population(pop_idx,29);
            
            % Calculate the new chisquared
            [chisquared,chisquared_array] = chiSqCalc(rates, A, P, sigma_A, C2_time, C4_time);
            %Add the chisquared value to its corresponding slot in the population
            pop_chisquared_array(pop_idx) = chisquared;
        end
            
        % Evaluate Generation
        % Sort from highest rated to lowest rated individuals
        [pop_chisquare_sorted,index_arr] = sort(pop_chisquared_array);
        population = population(index_arr,:);
        
        %Add to the best guess array the best guess of the population
        guess(genNum,1:Nparam) = population(1,:); % Current best guess
        guess(genNum,Nparam+1) = pop_chisquare_sorted(1); % Current best guess' fit value
        
        %population(1,:) % Display best guess to screen
        %     popfit(1);
        if guessUpdateMode == 1
            disp(['     Gen ' num2str(genNum) ': Best fitness guess so far: ' num2str(pop_chisquare_sorted(1))]);
        end
%% Part 4: genX
        %----------------------------------------------------------------------
        % (7) Plot the distribution of guesses (PART 4: genX)
        %----------------------------------------------------------------------
        if plot_GenerationalChiSquaredHistogramMode == 1
            figure(7)
            set(gcf,'Name','Chi-squared values of initial guesses');
            set(gcf,'Color','w');
            histogram(pop_chisquared_array);
            xlabel('\chi^2 Value');
            ylabel('Frequency');
            title(['Gen' num2str(genNum) ' Distribution of guesses chi-squared']);
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
        % * t32 will be determined by others in loop
        t24 = guess(genNum,6);
        t42 = guess(genNum,7);
        t25 = guess(genNum,8);
        t52 = guess(genNum,9);
        t56 = guess(genNum,10);
        t65 = guess(genNum,11);
        % * t63 will be determined by others in loop
        t36 = guess(genNum,12);
        t37 = guess(genNum,13);
        t73 = guess(genNum,14);
        t68 = guess(genNum,15);
        t86 = guess(genNum,16);
        
        % Loop condition
%         t32 = (t12 * t23 * t31)/(t13 * t21);
%         t63 = (t23 * t36 * t65* t52)/(t32 * t25 * t56);
        
        A1 = guess(genNum,17);        %21);
        A2 = guess(genNum,18);        %22);
        A3 = guess(genNum,19);        %23);
        A4 = guess(genNum,20);        %24);
        A5 = guess(genNum,21);        %25);
        A6 = guess(genNum,22);        %26);
        A7 = guess(genNum,23);        %27);
        A8 = guess(genNum,24);
        
        chisquared = guess(genNum,25);
        genNum_array = [genNum_array; genNum];
        chisquared_array = [chisquared_array; chisquared];
        
        % Redefine rates, kij, from the tij's
          k12 = 1/t12;
          k13 = 1/t13;
          k21 = 1/t21;
          k23 = 1/t23;
          k32 = 1/t32;
          k31 = 1/t31;
          k24 = 1/t24;
          k42 = 1/t42;
          k25 = 1/t25;
          k52 = 1/t52;
          k56 = 1/t56;
          k65 = 1/t65;
          k63 = 1/t63;
          k36 = 1/t36;
          k37 = 1/t37;
          k73 = 1/t73;
          k68 = 1/t68;
          k86 = 1/t86;
%           k89 = 1/t89;
%           k98 = 1/t98;
%           t96 = 1/t96;
%           t69 = 1/t69;

% Recalculate loop conditions
          k32 = (k12 * k23 * k31)/(k13 * k21);
          t32 = 1/k32;
          k63 = (k23 * k36 * k65* k52)/(k32 * k25 * k56);
          t63 = 1/k63;
          
        A = [A1;A2;A3;A4;A5;A6;A7;A8];
        rates = [k12,k13,k21,k23,k32,k31,k24,k42,k25,k52,k56,k65,k63,k36,k37,k73,k68,k86];

        if verboseMode == 1  
            fprintf(['Best fit from initial generation was member #%d:\n '...
                't12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
                't24 = %f, t42 = %f, t25 = %f, t52 = %f, t56 = %f, t65 = %f'...
                't63 = %f, t36 = %f,t37 = %f, t73 = %f,t68 = %f, t86 = %f'...
                '\n A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f,'...
                'A6 = %f, A7 = %f, A8 = %f,Best_chisquared = %f\r\n'],...
                index_arr(1),t12,t13,t21,t23,t31,t32,...
                t24, t42, t25, t52, t56, t65, t63, t36,...
                t37, t73, t68, t86,A1, A2, A3,A4, A5, A6, A7, A8,...
                Best_chisquared);
        end
        % Define K matrix
         % 8 state model:
         K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0;...
             k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0;...
             k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0;...
             0, k24, 0, -k42, 0, 0, 0, 0;...
             0, k25, 0, 0, -(k52+k56), k65, 0, 0;...
             0, 0, k36, 0, k56, -(k65+k63+k68), 0, k86;...
             0, 0, k37, 0, 0, 0, -k73, 0;...
             0, 0, 0, 0, 0, k68, 0, -k86;];
         
         %          % 9 state model:
         %          K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0, 0;...
         %              k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0, 0;...
         %              k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0, 0;...
         %              0, k24, 0, -k42, 0, 0, 0, 0, 0;...
         %              0, k25, 0, 0, -(k52+k56), k65, 0, 0, 0;...
         %              0, 0, k36, 0, k56, -(k65+k63+k68+k69), 0, k86, k96;...
         %              0, 0, k37, 0, 0, 0, -k73, 0, 0;...
         %              0, 0, 0, 0, 0, k68, 0, -(k86+k89),k98;...
         %              0, 0, 0, 0, 0, k69, 0, k89, -(k96+k98)];
         
         % Deine our new  P using this best guess K:
         [P, ~, ~, time] = k2P(K,time);
         
         %--------------------------------------------------------------------------
         % (1) Optimization Target #1: 1-D FRET HISTOGRAM (PART 4: genX)
         %--------------------------------------------------------------------------
         if fitHistMode == 1
              figure(1)
                hold on;
                [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate(P, A, sigma_A);
                
                if plotMode == 1
                            figure(1);
                            clf;
                            set(gcf,'Color','w');
                            set(gcf,'Name','FRET Histogram');
                            lineStyle = char('r--','g--','b--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                            
                            FRET_bins = linspace(0,1,100);
                            
                            if exist('data_hist_Plot','var') == 1
                                delete(data_hist_Plot)
%                                 disp('data_hist_Plot deleted')
                            end
                            % Replot data Histogram 
                            data_hist_Plot = plot(FRET_bins,targetHistogram);
                            data_hist_Plot.LineWidth = 2;
                            data_hist_Plot.Color = 'blue';
                            x.DisplayName = 'Data';
                            
                            xlabel('FRET Efficiency','FontSize',14);
                            ylabel('Frequency','FontSize',14);
                            
                            [sample_description, ~] = sample_descriptionGetter();
                            title_str = ['Experimental vs Simulated Histograms' ...
                                10 sample_description];
                            title(title_str,'fontsize',14);
                            
                            hold on;

                            % Clear plots from previous run
                            if exist('histPlot','var') == 1
                                delete(histPlot)
                                disp('histPlot deleted')
                            end
                            
                            % Plot histogram of each state
                            for i = 1:N
                                histPlot = plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
                                    lineStyle(i,:),'LineWidth',1,'DisplayName',Peq_names(i,:));
                                lgd = legend('show');
                                % lgd.Location = 'northwest';
                                lgd.FontSize = 14;
                                hold on
                            end
                            set(gca, 'FontSize', 14);
                            if exist('histPlotTot','var') == 1
                                delete(histPlotTot)
                            end
                            hold on
                            % Plot overall histogram
                            histPlotTot = plot(FRET_bins, hist_sim,...
                                lineStyle(i),'LineWidth',1,'DisplayName','Total Fit');
                            lgd_tot = legend('show');
                            
                end
         end
         %--------------------------------------------------------------------------
         % (2) Optimization Target #2: 2-point TCF (C2) (20 uSec and on) (PART 4: genX)
         %--------------------------------------------------------------------------
         if fitC2Mode == 1
             [time, C2, ~] = P2C(P, K, time);
                    if plotMode == 1
                        figure(2)
                        set(gcf,'Color','w');
                        
                        if exist('C2_plot','var')  == 1
                            delete(C2_plot)
                        end
                        
                        C2_plot = plot(time,C2);
                        title_str = ['Two point time correlation function'];
                        title(title_str,'FontSize',18);
                        xlabel('\tau (sec)','fontsize',16);
                        ylabel('C^{(2)}(\tau)','fontsize',16);
                        set(gca,'yscale','linear');
                        set(gca,'xscale','log');
                        set(gca,'FontSize',14);
                        grid on
                        axis tight;
                        
                        drawnow();
                    end
         end
         
         %--------------------------------------------------------------------------
         % (3) Optimization Target #3: 4-point TCF (C4) (10 uSec and on) (PART 4: genX)
         %--------------------------------------------------------------------------
         if fitC4Mode == 1
             [time,~,C4] = P2C(P, K, time);
             
             if plotMode == 1
                 figure(3)
                 
                 if exist('C4_plot','var')  == 1
                     delete(C4_plot)
                 end
                 
                 set(gcf,'Color','w');
                 hold on;
                 C4_plot = surf(time,time,C4);
                 title_str = ['Four point time correlation function'];
                 title(title_str,'FontSize',18);
                 xlabel('\tau_1 (sec)','fontsize',16);
                 ylabel('\tau_3 (sec)','fontsize',16);
                 zlabel('C^{(4)}(\tau)','fontsize',16);
                 set(gca,'yscale','log');
                 set(gca,'xscale','log');
                 set(gca,'FontSize',14);
                 grid on
                 axis tight;
                 
                 drawnow();
             end
         end
         %--------------------------------------------------------------------------
         % (4) Plotting the progress as a function of Generation
         %--------------------------------------------------------------------------
         figure(4);
         hold on;
         plot(genNum_array,chisquared_array,'o',...
             'MarkerSize',10,...
             'MarkerEdgeColor','b',...
             'MarkerFaceColor',[0.5,0.5,0.5]);
         drawnow();
         %------------------------------------------------------------------------------------------
         % (5) Update paramter-histograms as a function of the current best guesses  (PART 4: genX)
         %------------------------------------------------------------------------------------------
         if plot_paramater_histogram_mode == 1
             figure(5)
             param_strings = {'t12','t13','t21','t23','t31','t24','t42','t25','t52',...
                 't56','t65','t36','t37','t73','t68','t86','t32','t63',...
                 'A1','A2','A3','A4','A5','A6','A7','A8'};
             maxCountArray = zeros(1,Nparam);
             N_tij = Nparam - N; % N = number of states, N_tij is number of times *********************************
             for param_idx = 1:Nparam
                 subplot(2,4,param_idx)
                 
                 if ismember(param_idx,1:N_tij)                      %Rate histograms
                     nbins = 50;
                     histogram(population(:,param_idx),nbins);
                     [N,edges] = histcounts(population(:,param_idx),nbins);
                     maxCountArray(param_idx) = max(N);
                     xlabel('Time (Sec)');
                 elseif ismember(param_idx,(N_tij+1):Nparam) == 1    %FRET Historgrams
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
         end % if plot_parameter_histogram_mode
         
         %--------------------------------------------------------------------------
         % (6) Display a graphic of the network        (Part 5: Final Values)
         %--------------------------------------------------------------------------
         if plot_model_GenerationalUpdatesMode == 1
             
             figure(6);
             clf;
             set(gcf,'Name','Model: 8 state');
             
             xlim([0 4.1]);
             ylim([0 4.1]);
             
             %Designate spots for the states
             state1_loc = [0.5 2];
             state2_loc = [1.5 2.5];
             state3_loc = [1.5 1.5];
             state4_loc = [2.5 3.5];
             state5_loc = [2.5 2.5];
             state6_loc = [2.5 1.5];
             state7_loc = [2.5 0.5];
             state8_loc = [3.5 1.5];
             % or for 9 state
             %              state8_loc = [4 2.5];
             %              state9_loc = [4 1.5];
             hold on;
             text(0.1,4,['Gen# = ' num2str(genNum)],'FontSize',14);
             
            % Plot the state symbols
             text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
             text(state2_loc(1),state2_loc(2),'2','FontSize',24);
             text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
             text(state4_loc(1)-0.04,state4_loc(2)+0.05,'4','FontSize',24)
             text(state5_loc(1),state5_loc(2),'5','FontSize',24);
             text(state6_loc(1)-0.05,state6_loc(2),'6','FontSize',24);
             text(state7_loc(1)-0.04,state7_loc(2)+0.05,'7','FontSize',24)
             text(state8_loc(1),state8_loc(2),'8','FontSize',24);
%              text(state9_loc(1),state9_loc(2),'9','FontSize',24);

             
             %Plot the FRET values
             text(state1_loc(1),state1_loc(2)-0.2,['=' num2str(A1,'%.2f')],'FontSize',16, 'Color','b');
             text(state2_loc(1),state2_loc(2)-0.2,['=' num2str(A2,'%.2f')],'FontSize',16, 'Color','b');
             text(state3_loc(1),state3_loc(2)-0.2,['=' num2str(A3,'%.2f')],'FontSize',16, 'Color','b');
             text(state4_loc(1),state4_loc(2)-0.2,['=' num2str(A4,'%.2f')],'FontSize',16, 'Color','b');
             text(state5_loc(1),state5_loc(2)-0.2,['=' num2str(A5,'%.2f')],'FontSize',16, 'Color','b');
             text(state6_loc(1),state6_loc(2)-0.2,['=' num2str(A6,'%.2f')],'FontSize',16, 'Color','b');
             text(state7_loc(1),state7_loc(2)-0.2,['=' num2str(A7,'%.2f')],'FontSize',16, 'Color','b');
             text(state8_loc(1),state8_loc(2)-0.2,['=' num2str(A8,'%.2f')],'FontSize',16, 'Color','b');
%              text(state9_loc(1),state8_loc(2)-0.2,['=' num2str(A9,'%.2f')],'FontSize',16);

             
             %Plot Lines Between states
             line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
             line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
             line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
             line([state2_loc(1) state4_loc(1)],[state2_loc(2) state4_loc(2)],'Color','k','LineStyle','-');
             line([state2_loc(1) state5_loc(1)],[state2_loc(2) state5_loc(2)],'Color','k','LineStyle','-');
             line([state5_loc(1) state6_loc(1)],[state5_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
             line([state3_loc(1) state6_loc(1)],[state3_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
             line([state3_loc(1) state7_loc(1)],[state3_loc(2) state7_loc(2)],'Color','k','LineStyle','-');
             line([state6_loc(1) state8_loc(1)],[state6_loc(2) state8_loc(2)],'Color','k','LineStyle','-');
             
             %              % For 9 states
             %              line([state8_loc(1) state9_loc(1)],[state8_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
             %              line([state6_loc(1) state9_loc(1)],[state6_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
             %
             %Plot the inverse of the rates
             t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
             t12_msg = ['t_{12} = ' num2str(round(t12,6)) 10 't_{21} = ' num2str(round(t21,6))];
             text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',10);
             
             t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
             t23_msg = ['t_{23} = ' num2str(round(t23,6)) 10 't_{32} = ' num2str(round(t32,6))];
             text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',10);
             
             t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
             t31_msg = ['t_{31} = ' num2str(round(t31,6)) 10 't_{13} = ' num2str(round(t13,6))];
             text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',10);
             
             t24_loc = [state2_loc(1)+abs((state2_loc(1)- state4_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state4_loc(2))/2)];
             t24_msg = ['t_{24} = ' num2str(round(t24,6)) 10 't_{42} = ' num2str(round(t42,6))];
             text(t24_loc(1),t24_loc(2),t24_msg,'FontSize',10);
             
             t37_loc = [state7_loc(1) - abs((state3_loc(1) - state7_loc(1))/2),state7_loc(2)+ abs((state3_loc(2)- state7_loc(2))/2)];
             t37_msg = ['t_{37} = ' num2str(round(t37,6)) 10 't_{73} = ' num2str(round(t73,6))];
             text(t37_loc(1),t37_loc(2),t37_msg,'FontSize',10);
             
             t25_loc = [state2_loc(1)+abs((state2_loc(1)- state5_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state5_loc(2))/2)];
             t25_msg = ['t_{25} = ' num2str(round(t25,6)) 10 't_{52} = ' num2str(round(t52,6))];
             text(t25_loc(1),t25_loc(2),t25_msg,'FontSize',10);
             
             t36_loc = [state3_loc(1)+abs((state3_loc(1)- state6_loc(1))/2),state3_loc(2)+abs((state3_loc(2)- state6_loc(2))/2)];
             t36_msg = ['t_{36} = ' num2str(round(t36,6)) 10 't_{63} = ' num2str(round(t63,6))];
             text(t36_loc(1),t36_loc(2),t36_msg,'FontSize',10);
             
             t56_loc = [state5_loc(1)+abs((state5_loc(1)- state6_loc(1))/2),state5_loc(2)-abs((state5_loc(2)- state6_loc(2))/2)];
             t56_msg = ['t_{56} = ' num2str(round(t56,6)) 10 't_{65} = ' num2str(round(t65,6))];
             text(t56_loc(1),t56_loc(2),t56_msg,'FontSize',10);
             
             t68_loc = [state6_loc(1)+abs((state8_loc(1)- state6_loc(1))/2),state6_loc(2)+abs((state8_loc(2)- state6_loc(2))/2)];
             t68_msg = ['t_{68} = ' num2str(round(t36,6)) 10 't_{86} = ' num2str(round(t86,6))];
             text(t68_loc(1),t68_loc(2),t68_msg,'FontSize',10);
             
             %Plot the chi-squared value
             text(0.5,3.7,['\chi^2 = ' num2str(chisquared)],'FontSize',14);
             
             %Plot a title with the information in it
             [sample_description, save_prefix] = sample_descriptionGetter();
             title_str = ['8 state: ' sample_description];
             text(2.5,4,title_str,'fontsize',16);
             ylim([0,1.1]);
             
             axis off;
             
         end
        end
        %----------------------------------------------------------------------
        % Determine if the guesses have converged (PART 4: genX)
        %----------------------------------------------------------------------
        if (guess(genNum-1,Nparam+1) - guess(genNum,Nparam+1))/guess(genNum,Nparam+1) <= threshold
            Ntrials = Ntrials + 1;
        else
            Ntrials = 0;
        end
        if Ntrials == maxRepeats
            guess = guess(1:genNum,:);
            if verboseMode == 1
                disp(['The fitness has no longer improved after ' num2str(maxRepeats) ' repeats. Ending search.']);
            end
            break;
        end
        if clockMode == 1
            elapsedTime = toc;
            disp(['The amount of time per generation is: ' num2str(elapsedTime) ' seconds']);
        end
        if pauseBetweenGenerationMode == 1
            disp(['     Gen' num2str(genNum) ' Chi^2 = ' num2str(chisquared) ': Press Enter to proceed to the next generation']);
            pause();
        end
        end
%% Part 5: display results        
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
        % * t32 will be determined by others in loop
        t24 = guess(genNum,6);
        t42 = guess(genNum,7);
        t25 = guess(genNum,8);
        t52 = guess(genNum,9);
        t56 = guess(genNum,10);
        t65 = guess(genNum,11);
        % * t63 will be determined by others in loop
        t36 = guess(genNum,12);
        t37 = guess(genNum,13);
        t73 = guess(genNum,14);
        t68 = guess(genNum,15);
        t86 = guess(genNum,16);
        
        % Loop condition
%         t32 = (t12 * t23 * t31)/(t13 * t21);
%         t63 = (t23 * t36 * t65* t52)/(t32 * t25 * t56);
        
        A1 = guess(genNum,17);        %21);
        A2 = guess(genNum,18);        %22);
        A3 = guess(genNum,19);        %23);
        A4 = guess(genNum,20);        %24);
        A5 = guess(genNum,21);        %25);
        A6 = guess(genNum,22);        %26);
        A7 = guess(genNum,23);        %27);
        A8 = guess(genNum,24);
        
        chisquared = guess(genNum,25);
        genNum_array = [genNum_array; genNum];
        chisquared_array = [chisquared_array; chisquared];
        
        % Redefine rates, kij, from the tij's
          k12 = 1/t12;
          k13 = 1/t13;
          k21 = 1/t21;
          k23 = 1/t23;
          k32 = 1/t32;
          k31 = 1/t31;
          k24 = 1/t24;
          k42 = 1/t42;
          k25 = 1/t25;
          k52 = 1/t52;
          k56 = 1/t56;
          k65 = 1/t65;
          k63 = 1/t63;
          k36 = 1/t36;
          k37 = 1/t37;
          k73 = 1/t73;
          k68 = 1/t68;
          k86 = 1/t86;
%           k89 = 1/t89;
%           k98 = 1/t98;
%           t96 = 1/t96;
%           t69 = 1/t69;

% Recalculate loop conditions
          k32 = (k12 * k23 * k31)/(k13 * k21);
          t32 = 1/k32;
          k63 = (k23 * k36 * k65* k52)/(k32 * k25 * k56);
          t63 = 1/k63;
          
        A = [A1;A2;A3;A4;A5;A6;A7;A8];
        rates = [k12,k13,k21,k23,k32,k31,k24,k42,k25,k52,k56,k65,k63,k36,k37,k73,k68,k86];

        
        % Define K matrix
         % 8 state model:
         K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0;...
             k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0;...
             k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0;...
             0, k24, 0, -k42, 0, 0, 0, 0;...
             0, k25, 0, 0, -(k52+k56), k65, 0, 0;...
             0, 0, k36, 0, k56, -(k65+k63+k68), 0, k86;...
             0, 0, k37, 0, 0, 0, -k73, 0;...
             0, 0, 0, 0, 0, k68, 0, -k86;];
         
         %          % 9 state model:
         %          K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0, 0;...
         %              k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0, 0;...
         %              k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0, 0;...
         %              0, k24, 0, -k42, 0, 0, 0, 0, 0;...
         %              0, k25, 0, 0, -(k52+k56), k65, 0, 0, 0;...
         %              0, 0, k36, 0, k56, -(k65+k63+k68+k69), 0, k86, k96;...
         %              0, 0, k37, 0, 0, 0, -k73, 0, 0;...
         %              0, 0, 0, 0, 0, k68, 0, -(k86+k89),k98;...
         %              0, 0, 0, 0, 0, k69, 0, k89, -(k96+k98)];
         
         % Deine our new  P using this best guess K:
         [P, ~, ~, time] = k2P(K,time);    
        
         Peq = diag(P(:,:,end));
        if verboseMode == 1  
            fprintf(['Best fit from initial generation was member #%d:\n '...
                't12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
                't24 = %f, t42 = %f, t25 = %f, t52 = %f, t56 = %f, t65 = %f'...
                't63 = %f, t36 = %f,t37 = %f, t73 = %f,t68 = %f, t86 = %f'...
                '\n A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f,'...
                'A6 = %f, A7 = %f, A8 = %f,'...
                'P1eq = %f, P2eq = %f, P3eq = %f, P4eq = %f, P5eq = %f,'...
                'P6eq = %f, P7eq = %f, P8eq = %f,Best_chisquared = %f\r\n'],...
                index_arr(1),t12,t13,t21,t23,t31,t32,...
                t24, t42, t25, t52, t56, t65, t63, t36,...
                t37, t73, t68, t86,A1, A2, A3, A4, A5, A6, A7, A8,...
                Peq(1), Peq(2), Peq(3), Peq(4), Peq(5), Peq(6),...
                Peq(7), Peq(8),Best_chisquared);
        end
        
%% Part 5: Final values        
        %--------------------------------------------------------------------------
        % Save the data (PART 5: Final Values)
        %--------------------------------------------------------------------------
        
        %--------------------------------------------------------------------------
        % 5.1: Figure out where to save the data
        %--------------------------------------------------------------------------
        if saveMode == 1
            clear('outputFolderName')
            %Make output folder if it doesnt exist
            lowestChiSquareFolderName = [programName '_output' filesep() 'lowestChiSquare'];
            secondLowestChiSquareFolderName = [programName '_output' filesep() 'secondLowestChiSquare'];
            thirdLowestChiSquareFolderName = [programName '_output' filesep() 'thirdLowestChiSquare'];
            if exist(lowestChiSquareFolderName,'dir') ~= 7
                mkdir(lowestChiSquareFolderName);
                disp(['Making a folder called "' lowestChiSquareFolderName '" to hold the output']);
            end
            
            %If the folder exist and there are files in it, check the new fit
            %against them.
            bestfitFileName = [lowestChiSquareFolderName filesep() 'BestFitResults' '.mat'];
            if exist(bestfitFileName,'file') == 2
                save(bestfitFileName,'iter','-append');%***one place we update
                disp([save_prefix ': Updated the best fit with another iteration: #' num2str(iter)]);
                
                %Update the fitresults with the iteration number no matter what
                %Save the data in a higher chisquare folder
                outputFolderName = [programName '_output' filesep() 'higherChiSquares'];
                if exist(outputFolderName,'dir') ~= 7
                    mkdir(outputFolderName);
                    disp(['Making a folder called "' outputFolderName '" to hold the output']);
                end
                tempFitfoutName = ['fitResults_iter' num2str(iter,'%04d') '.mat'];
                disp(['Saving the results as ' tempFitfoutName '. ChiSqaure = ' num2str(chisquared)]);
                save([outputFolderName filesep() tempFitfoutName],...
                    'A1','A2','A3','A4','A5','A6','A7','A8','Peq',...
                    't12','t13','t21','t23','t31','t32',...
                    't24', 't42', 't25', 't52', 't56', 't65', 't63', 't36',...
                    't37', 't73', 't68', 't86',...
                    'chisquared','guess','genNum','iter',...
                    'k12','k13','k21','k23','k31','k32',...
                    'k24', 'k42', 'k25', 'k52', 'k56', 'k65', 'k63', 'k36',...
                    'k37', 'k73', 'k68', 'k86',...
                    'constructName','sample_description', 'save_prefix');
                
            
                current_LowestChisquared = chisquared;
                load(bestfitFileName,'chisquared','iter');
                prev_LowestChisquared = chisquared;%Save the previous X2 as prev_LowestChisquared
                chisquared = current_LowestChisquared;%reset chisquared to its current value
                if current_LowestChisquared < prev_LowestChisquared
                    %Just found a new best fit. Need to save it.
                    
                    fprintf('Congratulations! %f < %f (#1 fit). Will update best fit .\r\n',current_LowestChisquared,prev_LowestChisquared)
                    
                    %Moving the previous 2nd best fit to the number 3 slot
                    
                    if exist(secondLowestChiSquareFolderName,'dir') == 7
                        copyfile(secondLowestChiSquareFolderName,thirdLowestChiSquareFolderName);
                        disp('     Copied the previous second best fit into the third place slot');
                    end
                    
                    %Moving the previous best fit to the number 2 slot
                    copyfile(lowestChiSquareFolderName,secondLowestChiSquareFolderName);
                    disp('     Copied the previous best fit into the second place slot');
                    
                    outputFolderName = lowestChiSquareFolderName; %Assign the save folder to the lowest chi square value
                    
                elseif current_LowestChisquared == prev_LowestChisquared
                    
                    fprintf('Nothing changed.%f = %f. \r\n',current_LowestChisquared,prev_LowestChisquared);
                    saveMode = 0;
                    
                elseif current_LowestChisquared > prev_LowestChisquared
                    fprintf('Oh no! %f > %f (#1 fit). Will check to see if it is lower than #2 best fit .\n',current_LowestChisquared,prev_LowestChisquared)
                    % Two options:
                    % 1) The current x2 is higher than the number 2 and 3. Update iteration number on the best fit
                    % 2) current x2 is higher than #1 but less than #3. Save it as #2.
                    
                    %------------------------------------------------------------------
                    % Check if  the new chisquared deserves spot number 2
                    %------------------------------------------------------------------
                     if exist(secondLowestChiSquareFolderName,'dir') == 7
                        load([secondLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'chisquared');
                        save([secondLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'iter','-append');
                        prev_2ndBestChisquared = chisquared;
                        chisquared = current_LowestChisquared;%reset chisquared to its current value
                        if current_LowestChisquared < prev_2ndBestChisquared
                            fprintf('Congratulations! %f < %f (#2 fit). You found a new #2 fit.  Updating 2nd best fit.\r\n',current_LowestChisquared,prev_2ndBestChisquared)
                            %Move the number 2 fit to the number 3 slot
                            copyfile(secondLowestChiSquareFolderName,thirdLowestChiSquareFolderName);
                            disp('     Copied the previous second best fit into the third place slot');
                            outputFolderName = secondLowestChiSquareFolderName;
                        elseif current_LowestChisquared == prev_2ndBestChisquared
                            fprintf('Nothing changed.%f = %f (#2 fit).  \r\n',current_LowestChisquared,prev_2ndBestChisquared);
                            saveMode = 0;
                        elseif current_LowestChisquared > prev_2ndBestChisquared
                            fprintf('Oh No! %f > %f (#2 fit). Checking to see if fit is better than #3:\n',current_LowestChisquared,prev_2ndBestChisquared);
                            %------------------------------------------------------------------
                            % Check if  the new chisquared deserves spot number 3
                            %------------------------------------------------------------------
                            if exist(thirdLowestChiSquareFolderName,'dir') == 7
                                current_iter = iter;
                                load([thirdLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'chisquared','iter');
                                prev_iter = iter;
                                iter = current_iter;%Reset the iteration to the current one
                                save([thirdLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'iter','-append');
                                prev_3rdBestChisquared = chisquared;
                                chisquared = current_LowestChisquared;%reset chisquared to its current value;
                                if current_LowestChisquared < prev_3rdBestChisquared
                                    fprintf('Congratulations! %f < %f (#3fit). You found a new #3 fit.  Updating 3rd best fit.\r\n',current_LowestChisquared,prev_3rdBestChisquared)
                                    outputFolderName = thirdLowestChiSquareFolderName;
                                    copyfile([thirdLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],[programName '_output' filesep() 'higherChiSquares'  filesep() 'fitResults_iter' num2str(prev_iter,'%04d') '.mat']);
                                    
                                elseif current_LowestChisquared == prev_3rdBestChisquared
                                    fprintf('Nothing changed.%f = %f (#3 fit). \r\n',current_LowestChisquared,prev_3rdBestChisquared);
                                    saveMode = 0;
                                elseif current_LowestChisquared > prev_3rdBestChisquared
                                    fprintf('Oh No! %f > %f (#3 fit).\n',current_LowestChisquared,prev_3rdBestChisquared);
                                    
                                    saveMode = 0;
                                    
                                end
                            elseif exist(thirdLowestChiSquareFolderName,'dir') ~= 7
                                outputFolderName = thirdLowestChiSquareFolderName;
                                
                                mkdir(outputFolderName);
                                disp(['Making a folder called "' outputFolderName '" to hold the output']);
                                
                            end
                        end
                        elseif exist(secondLowestChiSquareFolderName,'dir') ~= 7
                        
                        outputFolderName = secondLowestChiSquareFolderName;
                        
                        mkdir(outputFolderName);
                        disp(['Making a folder called "' outputFolderName '" to hold the output']);
                        
                    end
                end
                elseif exist(foutName,'file') ~= 2
                iter = 1;
                outputFolderName = lowestChiSquareFolderName;
                
            end
            
        end
        
        %------------------------------------------------------------------
        % Update the folders with the correct information
        %------------------------------------------------------------------
        if saveMode == 1
            disp(['Will save the data in ' outputFolderName]);
            
            foutName = [outputFolderName filesep() 'fitInputData.mat'];
            if fitHistMode == 1
                save(foutName,'histogram_FilePath','FRET_bins','targetHistogram')
            end
            if fitC2Mode == 1
                if exist(foutName,'file') == 2
                    save(foutName,'C2_FilePath','C2_exp_x','C2_exp_y','yoff','weightC2func','-append');
                else
                    save(foutName,'C2_FilePath','C2_exp_x','C2_exp_y','yoff','weightC2func');
                end
                
            end
            if fitC4Mode == 1
                if exist(foutName,'file') == 2
                    save(foutName,'C4_FilePath','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff','-append');
                else
                    save(foutName,'C4_FilePath','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff');
                end
            end
            
            save([outputFolderName filesep() 'BestFitResults' '.mat'],...
                    'A1','A2','A3','A4','A5','A6','A7','A8','Peq',...
                    't12','t13','t21','t23','t31','t32',...
                    't24', 't42', 't25', 't52', 't56', 't65', 't63', 't36',...
                    't37', 't73', 't68', 't86',...
                    'k12','k13','k21','k23','k31','k32',...
                    'k24', 'k42', 'k25', 'k52', 'k56', 'k65', 'k63', 'k36',...
                    'k37', 'k73', 'k68', 'k86',...
                    'chisquared','guess','genNum','iter',...
                    'constructName','sample_description', 'save_prefix','yoff');
            %Save the paramaters used for the genetic algorithm
            save([outputFolderName filesep()  'genAlgParamaters.mat'],'programName','threshold',...
                'NmembersInitPop','maxGenerations','maxRepeats','Nreproduce','maxIterations',...
                'Nmutations','boundsArray','maxMutationArray',...
                'fitC2Mode','fitC4Mode','fitHistMode','fitHistData_mode',...
                'constructName','sample_description', 'save_prefix');
            %Save paramatres used to plot
            save([outputFolderName filesep() 'plottingParamaters.mat'],'FRET_bins','sigma_A',...
                'C2_exp_x','C2_exp_y','constructName','sample_description', 'save_prefix');
            
        end
        
        %--------------------------------------------------------------------------
        % Save the data PLOTS     (PART 5: Final Values)
        %--------------------------------------------------------------------------
        if plotMode == 1
            if showBestFitMode == 1
                disp('     Loading the best fit from all runs.');
                load([lowestChiSquareFolderName filesep() 'BestFitResults.mat'],...
                    'A1','A2','A3','A4','A5','A6','A7','A8','Peq',...
                    't12','t13','t21','t23','t31','t32',...
                    't24', 't42', 't25', 't52', 't56', 't65', 't63', 't36',...
                    't37', 't73', 't68', 't86',...
                    'k12','k13','k21','k23','k31','k32',...
                    'k24', 'k42', 'k25', 'k52', 'k56', 'k65', 'k63', 'k36',...
                    'k37', 'k73', 'k68', 'k86',...
                    'chisquared','guess','iter');
                saveMode = 1;
                outputFolderName = lowestChiSquareFolderName;
            end
            %//////////////////////////////////////////////////////////////////////////
            % Plot final results of the algorithm   (Part 5: Final Values)
            %//////////////////////////////////////////////////////////////////////////
            
            
            % Define K matrix
            % 8 state model:
            K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0;...
                k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0;...
                k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0;...
                0, k24, 0, -k42, 0, 0, 0, 0;...
                0, k25, 0, 0, -(k52+k56), k65, 0, 0;...
                0, 0, k36, 0, k56, -(k65+k63+k68), 0, k86;...
                0, 0, k37, 0, 0, 0, -k73, 0;...
                0, 0, 0, 0, 0, k68, 0, -k86;];
            
            %          % 9 state model:
            %          K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0, 0;...
            %              k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0, 0;...
            %              k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0, 0;...
            %              0, k24, 0, -k42, 0, 0, 0, 0, 0;...
            %              0, k25, 0, 0, -(k52+k56), k65, 0, 0, 0;...
            %              0, 0, k36, 0, k56, -(k65+k63+k68+k69), 0, k86, k96;...
            %              0, 0, k37, 0, 0, 0, -k73, 0, 0;...
            %              0, 0, 0, 0, 0, k68, 0, -(k86+k89),k98;...
            %              0, 0, 0, 0, 0, k69, 0, k89, -(k96+k98)];
            
            % Deine our new  P using this best guess K:
            [P, ~, ~, time] = k2P(K,time);
            
            Peq = diag(P(:,:,end));
            %---------------------------------------------------------------------------------------
            % (1) Optimization Target #1: Histogram  (Part 5: Final Values)
            %----------------------------------------------------------------------------------------
            %     if fitHistMode == 1
            [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate(P, A, sigma_A);
            
            if plotMode == 1
                figure(1);
                clf;
                set(gcf,'Color','w');
                set(gcf,'Name','FRET Histogram');
                lineStyle = char('r--','g--','b--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                
                FRET_bins = linspace(0,1,100);
                
                if exist('data_hist_Plot','var') == 1
                    delete(data_hist_Plot)
                    %                                 disp('data_hist_Plot deleted')
                end
                % Replot data Histogram
                data_hist_Plot = plot(FRET_bins,targetHistogram);
                data_hist_Plot.LineWidth = 2;
                data_hist_Plot.Color = 'blue';
%                 x.DisplayName = 'Data';
                            
                            xlabel('FRET Efficiency','FontSize',14);
                            ylabel('Frequency','FontSize',14);
                            
                            [sample_description, ~] = sample_descriptionGetter();
                            title_str = ['Experimental vs Simulated Histograms' ...
                                10 sample_description];
                            title(title_str,'fontsize',14);
                            
                            hold on;

                            % Clear plots from previous run
                            if exist('histPlot','var') == 1
                                delete(histPlot)
                                disp('histPlot deleted')
                            end
                            
                            % Plot histogram of each state
                            for i = 1:N
                                histPlot = plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
                                    lineStyle(i,:),'LineWidth',1,'DisplayName',Peq_names(i,:));
                                lgd = legend('show');
                                % lgd.Location = 'northwest';
                                lgd.FontSize = 14;
                                hold on
                            end
                            set(gca, 'FontSize', 14);
                            if exist('histPlotTot','var') == 1
                                delete(histPlotTot)
                            end
                            hold on
                            % Plot overall histogram
                            histPlotTot = plot(FRET_bins, hist_sim,...
                                lineStyle(i),'LineWidth',1,'DisplayName','Total Fit');
                            lgd_tot = legend('show');
            end  
                            
                            %---------------------------------------------------------------------------------------
                            % (2) Optimization Target #2: 2-point TCF (C2)  (20 uSec and on)  (Part 5: Final Values)
                            %----------------------------------------------------------------------------------------
%                             if fitC2Mode == 1
                                [time, C2, ~] = P2C(P, K, time);
                                if plotMode == 1
                                    figure(2)
                                    set(gcf,'Color','w');
                                    
                                    if exist('C2_plot','var')  == 1
                                        delete(C2_plot)
                                    end
                                    
                                    C2_plot = plot(time,C2);
                                    title_str = ['Two point time correlation function'];
                                    title(title_str,'FontSize',18);
                                    xlabel('\tau (sec)','fontsize',16);
                                    ylabel('C^{(2)}(\tau)','fontsize',16);
                                    set(gca,'yscale','linear');
                                    set(gca,'xscale','log');
                                    set(gca,'FontSize',14);
                                    grid on
                                    axis tight;
                                    
                                    drawnow();
                                end
                                % end
                                %---------------------------------------------------------------------------------------
                                % (3) Optimization Target #3: 4-point TCF (C4) (20 uSec and on)  (Part 5: Final Values)
                                %---------------------------------------------------------------------------------------
%                                 if fitC4Mode == 1
                                    [time,~,C4] = P2C(P, K, time);
                                    
                                    if plotMode == 1
                                        figure(3)
                                        
                                        if exist('C4_plot','var')  == 1
                                            delete(C4_plot)
                                        end
                                        
                                        set(gcf,'Color','w');
                                        hold on;
                                        C4_plot = surf(time,time,C4);
                                        title_str = ['Four point time correlation function'];
                                        title(title_str,'FontSize',18);
                                        xlabel('\tau_1 (sec)','fontsize',16);
                                        ylabel('\tau_3 (sec)','fontsize',16);
                                        zlabel('C^{(4)}(\tau)','fontsize',16);
                                        set(gca,'yscale','log');
                                        set(gca,'xscale','log');
                                        set(gca,'FontSize',14);
                                        grid on
                                        axis tight;
                                        
                                        drawnow();
                                    end
                                    %                                 end
                                    %--------------------------------------------------------------------------
                                    % (6) Display a graphic of the network        (Part 5: Final Values)
                                    %--------------------------------------------------------------------------
                                    figure(6);
                                    clf;
                                    set(gcf,'Name','Model: 8 state');
                                    
                                    xlim([0 4.1]);
                                    ylim([0 4.1]);
                                    
                                    %Designate spots for the states
                                    state1_loc = [0.5 2];
                                    state2_loc = [1.5 2.5];
                                    state3_loc = [1.5 1.5];
                                    state4_loc = [2.5 3.5];
                                    state5_loc = [2.5 2.5];
                                    state6_loc = [2.5 1.5];
                                    state7_loc = [2.5 0.5];
                                    state8_loc = [3.5 1.5];
                                    % or for 9 state
                                    %              state8_loc = [4 2.5];
                                    %              state9_loc = [4 1.5];
                                    hold on;
                                    text(0.1,4,['Gen# = ' num2str(genNum)],'FontSize',14);
                                    
                                    % Plot the state symbols
                                    text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
                                    text(state2_loc(1),state2_loc(2),'2','FontSize',24);
                                    text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
                                    text(state4_loc(1)-0.04,state4_loc(2)+0.05,'4','FontSize',24)
                                    text(state5_loc(1),state5_loc(2),'5','FontSize',24);
                                    text(state6_loc(1)-0.05,state6_loc(2),'6','FontSize',24);
                                    text(state7_loc(1)-0.04,state7_loc(2)+0.05,'7','FontSize',24)
                                    text(state8_loc(1),state8_loc(2),'8','FontSize',24);
                                    %              text(state9_loc(1),state9_loc(2),'9','FontSize',24);
                                    
                                    
                                    %Plot the FRET values
                                    text(state1_loc(1),state1_loc(2)-0.2,['=' num2str(A1,'%.2f')],'FontSize',16, 'Color','b');
                                    text(state2_loc(1),state2_loc(2)-0.2,['=' num2str(A2,'%.2f')],'FontSize',16, 'Color','b');
                                    text(state3_loc(1),state3_loc(2)-0.2,['=' num2str(A3,'%.2f')],'FontSize',16, 'Color','b');
                                    text(state4_loc(1),state4_loc(2)-0.2,['=' num2str(A4,'%.2f')],'FontSize',16, 'Color','b');
                                    text(state5_loc(1),state5_loc(2)-0.2,['=' num2str(A5,'%.2f')],'FontSize',16, 'Color','b');
                                    text(state6_loc(1),state6_loc(2)-0.2,['=' num2str(A6,'%.2f')],'FontSize',16, 'Color','b');
                                    text(state7_loc(1),state7_loc(2)-0.2,['=' num2str(A7,'%.2f')],'FontSize',16, 'Color','b');
                                    text(state8_loc(1),state8_loc(2)-0.2,['=' num2str(A8,'%.2f')],'FontSize',16, 'Color','b');
                                    %              text(state9_loc(1),state8_loc(2)-0.2,['=' num2str(A9,'%.2f')],'FontSize',16);
                                    
                                    
                                    %Plot Lines Between states
                                    line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
                                    line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
                                    line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
                                    line([state2_loc(1) state4_loc(1)],[state2_loc(2) state4_loc(2)],'Color','k','LineStyle','-');
                                    line([state2_loc(1) state5_loc(1)],[state2_loc(2) state5_loc(2)],'Color','k','LineStyle','-');
                                    line([state5_loc(1) state6_loc(1)],[state5_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
                                    line([state3_loc(1) state6_loc(1)],[state3_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
                                    line([state3_loc(1) state7_loc(1)],[state3_loc(2) state7_loc(2)],'Color','k','LineStyle','-');
                                    line([state6_loc(1) state8_loc(1)],[state6_loc(2) state8_loc(2)],'Color','k','LineStyle','-');
                                    
                                    %              % For 9 states
                                    %              line([state8_loc(1) state9_loc(1)],[state8_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
                                    %              line([state6_loc(1) state9_loc(1)],[state6_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
                                    %
                                    %Plot the inverse of the rates
                                    t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
                                    t12_msg = ['t_{12} = ' num2str(round(t12,6)) 10 't_{21} = ' num2str(round(t21,6))];
                                    text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',10);
                                    
                                    t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
                                    t23_msg = ['t_{23} = ' num2str(round(t23,6)) 10 't_{32} = ' num2str(round(t32,6))];
                                    text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',10);
                                    
                                    t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
                                    t31_msg = ['t_{31} = ' num2str(round(t31,6)) 10 't_{13} = ' num2str(round(t13,6))];
                                    text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',10);
                                    
                                    t24_loc = [state2_loc(1)+abs((state2_loc(1)- state4_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state4_loc(2))/2)];
                                    t24_msg = ['t_{24} = ' num2str(round(t24,6)) 10 't_{42} = ' num2str(round(t42,6))];
                                    text(t24_loc(1),t24_loc(2),t24_msg,'FontSize',10);
                                    
                                    t37_loc = [state7_loc(1) - abs((state3_loc(1) - state7_loc(1))/2),state7_loc(2)+ abs((state3_loc(2)- state7_loc(2))/2)];
                                    t37_msg = ['t_{37} = ' num2str(round(t37,6)) 10 't_{73} = ' num2str(round(t73,6))];
                                    text(t37_loc(1),t37_loc(2),t37_msg,'FontSize',10);
                                    
                                    t25_loc = [state2_loc(1)+abs((state2_loc(1)- state5_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state5_loc(2))/2)];
                                    t25_msg = ['t_{25} = ' num2str(round(t25,6)) 10 't_{52} = ' num2str(round(t52,6))];
                                    text(t25_loc(1),t25_loc(2),t25_msg,'FontSize',10);
                                    
                                    t36_loc = [state3_loc(1)+abs((state3_loc(1)- state6_loc(1))/2),state3_loc(2)+abs((state3_loc(2)- state6_loc(2))/2)];
                                    t36_msg = ['t_{36} = ' num2str(round(t36,6)) 10 't_{63} = ' num2str(round(t63,6))];
                                    text(t36_loc(1),t36_loc(2),t36_msg,'FontSize',10);
                                    
                                    t56_loc = [state5_loc(1)+abs((state5_loc(1)- state6_loc(1))/2),state5_loc(2)-abs((state5_loc(2)- state6_loc(2))/2)];
                                    t56_msg = ['t_{56} = ' num2str(round(t56,6)) 10 't_{65} = ' num2str(round(t65,6))];
                                    text(t56_loc(1),t56_loc(2),t56_msg,'FontSize',10);
                                    
                                    t68_loc = [state6_loc(1)+abs((state8_loc(1)- state6_loc(1))/2),state6_loc(2)+abs((state8_loc(2)- state6_loc(2))/2)];
                                    t68_msg = ['t_{68} = ' num2str(round(t36,6)) 10 't_{86} = ' num2str(round(t86,6))];
                                    text(t68_loc(1),t68_loc(2),t68_msg,'FontSize',10);
                                    
                                    %Plot the chi-squared value
                                    text(0.5,3.7,['\chi^2 = ' num2str(chisquared)],'FontSize',14);
                                    
                                    %Plot a title with the information in it
                                    [sample_description, save_prefix] = sample_descriptionGetter();
                                    title_str = ['8 state: ' sample_description];
                                    text(2.5,4,title_str,'fontsize',16);
                                    ylim([0,1.1]);
                                    
                                    axis off;
                                    
                                    
                                    
        end
    end
    iter = iter + 1;
    %--------------------------------------------------------------------------
    %The end of the program
    %--------------------------------------------------------------------------
    disp(['iter = ' num2str(iter) ': ' programName ' complete.']);

end



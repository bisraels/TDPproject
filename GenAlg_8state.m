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

% constructFolderNames = {'S4S5'};
constructFolderNames = {'S1S2'};
% constructFolderNames = {'S18S2'};
% constructFolderNames =  {'S19S5'};
% constructFolderNames = {'S19S5';'S1S2';'S4S5';'S18S2';};
%--------------------------------------------------------------------------
% Genetic Algorithm Patamaters
%--------------------------------------------------------------------------
NmembersInitPop = 100;          % Size of inital pop
maxGenerations = 1000;
maxRepeats = 10;
percentReproduce = 25;          % top 25% combine together (can be mutated)
percentToKeep = 50;             % Number of population to pass to next gen (can be mutated)
Nmutations = 1*NmemberInitPop;  % Number of mutations per generation
threshold = 0.001;              % Minimal allowable percentage  difference before program quits.
maxIterations = 1006;           % How is this different than maxGenerations?
forceMoreIterationsMode = 1;    % Independent iterations
NumberToExtend = 1;             
% NOTES: lowering the percent that reproduces but raising the percent to
% keep as high seems to help.
    % Why? What are these really doing in the GenAlg process?
    
%-------------------DERIVED PARAMATERS------------------
Nreproduct = round(NmembersInitPop * percentReproduce/100);  % # individuals who reproduce per generation.
membersToKeep = NmembersInitPop * percentToKeep/100;

%--------------------------------------------------------------------------
% Declare global variables
%--------------------------------------------------------------------------
global normalizeMode verboseMode guessUpdateMode diagnoseMode useAnalyticalAlgorithmsMode plotMode
global genNum fitHistMode fitC2Mode fitC4Mode
global targetHistogram weightingFactor_FREThist FRET_bins sigma_A1 sigma_A2 sigma_A3 sigma_A4 sigma_A5 sigma_A6 sigma_A7 sigma_A8 sigma_A9  % Histogram optimization
global C2_exp_x C2_exp_y weightingFactor_C2  weightC2func yoff
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func zoff

%--------------------------------------------------------------------------
% User Options                                                  (Part 1)
%--------------------------------------------------------------------------
normalizeMode = 0;      % Normalizes the 2pt and 4pt correlation functions
verboseMode = 1;
guessUpdateMode = 0;    % Very detailed: shows changes to any paramters
diagnoseMode = 0;       % Shows relative contribution of the 2pt TCF to the histogram at the end of each chi-squared calculation

clockMode = 0;          % Times various features of code.
saveModeDefault = 1;
plotMode = 1;           % Makes plots.
useFigurePosnMode = 1;

plot_Gen1_memberGuessesMode = 0; % Figure(1-2-3): Plots Histograms, C2, C4 of Gen1 (SLOW!)
plot_Gen1ChiSquaredMode = 0;     % Figure(4): Adds a point to the chi-sq vs member plot (plot4) (SLOW); Plot statistics from the first generation (informative)

plot_GenerationalProgressMode = 1; % Figure 4: Gray circles on Fig(4); Plot a generational update with best fit from the current generation
plot_parameter_histogram_mode = 0; % Figure 5

plot_model_GenerationalUpdatesMode = 1;       % Figure 6: Makes a model of the system with rates between states EACH GENERATIOON
plot_GenerationalChiSquaredHistogramMode = 0; % Figure 7: Histogram of chi-squared values
pauseBetweenGenerationMode = 0;               % Requires manual button press to continue

% At the end of the program
showBestFitMode = 1;  % Replaces current guess with best guess (end of code)

fitHistMode = 1;  fitHistData_mode = 1; % If 0 you will fit to the histogram fit
fitC2Mode = 1;
fitC4Mode = 1;

% Choose Weighting Amount
weightingFactor_FREThist = 50;  % Weighting for FRET hist comparison
weightingFactor_C2 = 1;
weightingFactor_C4_t0 = 1000;

if sum([fitHistMode,fitC2Mode,fitC4Mode]) == 0
    error('No surfaces to optimize to! Pick at least one.')
end

%**************************************************************************
useAnalyticalAlgorithmsMode = 1;
if useAnalyticalAlgorithmsMode == 1
    if verboseMode == 1
        disp('     ***using ANALYTICAL algorithms');
    end
else
    if verboseMode == 1
        disp('     ***using NUMERICAL algorithms');
        disp('           ==> Loading the conditional probabilities as a function of rates');   % Change this to a function of eigenvector components?
    end
    load('filename')
end

%**************************************************************************
ContributionArray = zeros(1,3);

if verboseMode == 1
    fprintf('       Each generation will have %d/%d of its members reproduce(%d%%)\r',Nreproduce,NmembersInitPop, percentReproduce);
    fprintf('       Each generation will see a total of %d mutations. (%d per individual on average).\r', Nmutations,Nmutations/NmembersInitPop);
    fprintf('       Each generation will have the top %d percent unmutated (%d/%d individuals).', percentToKeep,membersToKeep, NmembersInitPop);
end

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
t89_bounds = [1e-6, 1e-1];
t98_bounds = [1e-6, 1e-1];
% * t96_bounds will be determined by others in loop
t69_bounds = [1e-6, 1e-1];

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
    t24_bounds; t42_bounds; t25_bounds; t52_bounds; t56_bounds t63_bounds;...
    t36_bounds; t37_bounds; t73_bounds; t68_bounds; t86_bounds; t89_bounds;...
    t98_bounds; t69_bounds;...
    A1_bounds; A2_bounds; A3_bounds; A4_bounds; A5_bounds;A6_bounds; A7_bounds;...
    A8_bounds; A9_bounds];

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
t63_mutate = Max_mut_factor * t67_bounds(2)/2;
t36_mutate = Max_mut_factor * t76_bounds(2)/2;
t37_mutate = Max_mut_factor * t37_bounds(2)/2; 
t73_mutate = Max_mut_factor * t73_bounds(2)/2;
t68_mutate = Max_mut_factor * t68_bounds(2)/2;
t86_mutate = Max_mut_factor * t86_bounds(2)/2; 
t89_mutate = Max_mut_factor * t89_bounds(2)/2;
t98_mutate = Max_mut_factor * t98_bounds(2)/2;
% * t96 -> detailed balance
t69_mutate = Max_mut_factor * t69_bounds(2)/2; 

A1_mutate = (A1_bounds(2) - A1_bounds(1))*0.1;%0.05;
A2_mutate = (A2_bounds(2) - A2_bounds(1))*0.1;%0.05;0.05;
A3_mutate = (A3_bounds(2) - A3_bounds(1))*0.1;%0.05;0.05;
A4_mutate = (A4_bounds(2) - A4_bounds(1))*0.1;
A5_mutate = (A5_bounds(2) - A5_bounds(1))*0.1;
A6_mutate = (A6_bounds(2) - A6_bounds(1))*0.1;
A7_mutate = (A7_bounds(2) - A7_bounds(1))*0.1;
A8_mutate = (A8_bounds(2) - A8_bounds(1))*0.1;
A9_mutate = (A9_bounds(2) - A9_bounds(1))*0.1;

maxMutationArray = [t12_mutate;t13_mutate;t21_mutate;t23_mutate;t31_mutate;...
    t24_mutate; t42_mutate; t25_mutate; t52_mutate; t56_mutate; ...
    t63_mutate; t36_mutate; t37_mutate; t73_mutate; ...
    t68_mutate; t86_mutate; t89_mutate; t98_mutate; t69_mutate;...
    A1_mutate;A2_mutate;A3_mutate; A4_mutate; A5_mutate; A6_mutate; ...
    A7_mutate; A8_mutate; A9_mutate];

%Start in the single molecule folder (smData_Processed): comp specific
[computer_terminal_str, terminalID] = computerMode(pwd);

NconstructFolderNames = numel(constructFolderNames);

for construct_idx = 1:NconstructFolderNames
    constructName = char(constructFolderNames(construct_idx));
    disp(['     Construct' num2str(construct_idx) '/' num2str(NconstructFolderNames) ': ' constructName];
    
    if exist('constructName', 'var') ~= 1
        error('Not sure what construct to analyze');
    end
    
    % Find GenAlg files and load target data:
    [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix] = fileFinder(terminalID, constructName, modelName, protein_str);

    % Load file path for data plots
    [histogram_FilePath, C2_FilePath, C4_FilePath] = loadFilePath_dataPlots(genAlgDataLocationPrefix, protein_str);   
  
  %--------------------------------------------------------------------------
    % Set up the figure positions once and done
    %--------------------------------------------------------------------------
    if plotMode == 1
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
        data_hist_Plot.DisplayName = 'Data';
        
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
    end
    % end
    
    %--------------------------------------------------------------------------
    % (2) Optimization Target #2: 2-point TCF (C2)           (20 uSec and on)  (PART 2: Load the target data)
    %--------------------------------------------------------------------------
    if fitC2Mode == 1
        load(C2_FilePath,'time','yData','y','yoff');
        fitC2Data_mode = 1; % If 0 you will fit to the histogram fit
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
        
        load(C4_FilePath,'C4','tau1arraySec','tau3arraySec','tau2ValSec');
        C4_tau2eq0_exp = C4;
        C4_tau1range = tau1arraySec;
        C4_tau3range = tau3arraySec';
        tau2 = tau2ValSec;
        
        load(C2_FilePath,'time','yData','y','yoff');
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
        for pop_idx = 1:NmembersInitPop
            t12 = population(pop_idx,1);  k12 = 1/t12;
            t13 = population(pop_idx,2);
            t21 = population(pop_idx,3);
            t23 = population(pop_idx,4);
            t31 = population(pop_idx,5);
            t24 = population(pop_idx,6);
            t42 = population(pop_idx,7);
            t25 = population(pop_idx,8);
            t52 = population(pop_idx,9);
            t56 = population(pop_idx,10);
            % * t65_bounds will be determined by others in loop
            t63 = population(pop_idx,11);
            t36 = population(pop_idx,12);
            t37 = population(pop_idx,13);
            t73 = population(pop_idx,14);
            t68 = population(pop_idx,15);
            t86 = population(pop_idx,16);
%             t89 = population(pop_idx,17);
%             t98 = population(pop_idx,18);
            % * t96_bounds will be determined by others in loop
%             t69 = population(pop_idx,20);
            
            A1 = population(pop_idx,21);
            A2 = population(pop_idx,22);
            A3 = population(pop_idx,23);
            A4 = population(pop_idx,24);
            A5 = population(pop_idx,25);
            A6 = population(pop_idx,26);
            A7 = population(pop_idx,27);
            A8 = population(pop_idx,28);
%             A9 = population(pop_idx,29);
            
          if plotMode == 1
                %--------------------------------------------------------------------------
                % gen1: Plot the array of initial guesses             (PART 3: Gen1)
                %--------------------------------------------------------------------------
                if plot_Gen1_memberGuessesMode == 1
                    %--------------------------------------------------------------------------
                    % (1) Optimization Target #1: 1-D FRET HISTOGRAM  (PART 3: Gen1)
                    %--------------------------------------------------------------------------
                    if fitHistMode == 1
                        Peq = histMaker_Nstate(p, A, sigma_A);




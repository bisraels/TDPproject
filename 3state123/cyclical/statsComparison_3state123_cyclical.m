% function statsComparison_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3)
programName = 'statsComparison_3state123_cyclical';
close all;
%--------------------------------------------------------------------------
%USER OPTIONS
%--------------------------------------------------------------------------
saveMode = 0;
clockMode = 1;%if clockMode == 1, tic; end if clockMode == 1, disp(['     Took ' num2str(toc) ' seconds to run ' display_str]); end
plotMode = 1;

compareHistMode = 1;
compareC2Mode = 0;
compareC4Mode = 0;

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


%//////////////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------
% SET PARAMATERS
%--------------------------------------------------------------------------
% switch nargin
%     case 0
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
% end

% For testing (Has sensitive rates that will make a wrong solution diverge)
% t12 = 0.000825; t13 = 0.006767; t21 = 0.000615; t23 = 0.006725; t31 = 0.003205; t32 = 0.004270; A1 = 0.737326; A2 = 0.518160; A3 = 0.315085;

%--------------------------------------------------------------------------
% Set the rates (5 rates)
%--------------------------------------------------------------------------
k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;

k32 = k12*k23*k31/(k13*k21); % Detailed balance condition: %k31 will be the rate fixed by the others
t32 = 1/k32;

%//////////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------
% Compare the histograms
%--------------------------------------------------------------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if compareHistMode == 1
    disp('Comparing the time it takes to calculate histograms...');
    figure(21);
    clf
    set(gcf,'Color','w');
    
    sigma_A1 = 0.15;
    sigma_A2 = 0.1;
    sigma_A3 = 0.1;
    FRET_bins = linspace(0,1,100);
    
    
    
    %----------------------------------------------------------------------
    % Calculate the histogram using the analytical code             (Part3)
    %----------------------------------------------------------------------
    if clockMode == 1, tic; end
    [Peq] = histMaker_3state123_cyclical_analytical(t12,t13,t21,t23,t31,A1,A2,A3);
    display_str = 'histMaker_3state123_cyclical_analytical';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    p1_eq = Peq(1);
    p2_eq = Peq(2);
    p3_eq = Peq(3);
    hist_sim3 = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim3 = sum(hist_sim3);
    hist_sim3 = hist_sim3./denom_hist_sim3;
    fitHistPlot2 = plot(FRET_bins,hist_sim3);
    fitHistPlot2.LineStyle = ':';
    fitHistPlot2.Color = 'g';
    fitHistPlot2.LineWidth = 2;
    fitHistPlot2.DisplayName = ['histMaker 3state123 cyclical analytical'];
    hold on;
    
    %----------------------------------------------------------------------
    % Calculate the histogram using the new code
    %----------------------------------------------------------------------
    hold on
    if clockMode == 1
        tic;
    end
    [Peq] = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
    display_str = 'histMaker_3state123_cyclical';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    p1_eq = Peq(1);
    p2_eq = Peq(2);
    p3_eq = Peq(3);
    hist_sim2 = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim2 = sum(hist_sim2);
    hist_sim2 = hist_sim2./denom_hist_sim2;
    fitHistPlot2 = plot(FRET_bins,hist_sim2);
    fitHistPlot2.LineStyle = '--';
    fitHistPlot2.Color = 'blue';
    fitHistPlot2.LineWidth = 2;
    fitHistPlot2.DisplayName = ['histMaker 3state123 cyclical'];
    
    leg = legend;
    leg.Location = 'best';
    
    
    %--------------------------------------------------------------------------
    % Save the data
    %--------------------------------------------------------------------------
    if saveMode == 1
        %--------------------------------------------------------------------------
        % Make a folder to hold output
        %-------------------------------------------------------------------------
        %Make output folder if it doesnt exist
        outputFolderName = [programName '_output'];
        if exist(outputFolderName,'dir') ~= 7
            mkdir(outputFolderName);
            disp('Making a folder to hold the output');
        end
        
        %--------------------------------------------------------------------------
        % Save a figure
        %-------------------------------------------------------------------------
        set(gcf,'PaperPositionMode','auto');
        extension = 'png';
        runNum = 1;
        foutName = ['histSimulationComparison_' num2str(runNum)];
        filePathOut = [outputFolderName filesep() foutName];
        while exist([filePathOut '.' extension],'file') == 2
            runNum = runNum + 1;
            fprintf('Increasing the run num to %f',runNum);
            foutName = ['histSimulationComparison_' num2str(runNum)];
            filePathOut = [outputFolderName filesep() foutName];
        end
        print(gcf,filePathOut,['-d' extension],'-r0');
        disp(['Saving the figure as ' foutName]);
    end
end



%//////////////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------
% Compare the 2-point TCFs
%--------------------------------------------------------------------------
%//////////////////////////////////////////////////////////////////////////
if compareC2Mode == 1
    disp('Comparing the time it takes to calculate 2-point TCFs...');
    figure(22)
    clf
    set(gcf,'Color','w');
    
    %Make an array of time for the x-axis
    Npts = 150;
    time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
    C2_exp_x = time;
    
    
    
    %----------------------------------------------------------------------
    %Calculate the TCF with the old analytical expressions (cleaner)
    %----------------------------------------------------------------------
    if clockMode == 1
        tic;
    end
    %  tcf = C2maker_3state123_cyclical_analytical(t12,t13,t21,t23,t31,A1,A2,A3,time)
    tcf_sim2 = C2maker_3state123_cyclical_analytical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
    display_str = 'C2maker 3state123 cyclical analytical';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    tcf_plot = plot(C2_exp_x,tcf_sim2,'g:','LineWidth',2,'DisplayName',display_str);
    hold on;
    %RESULT:      Took 0.0043612 seconds to run C2maker 3state123 cyclical analytical
    
    %----------------------------------------------------------------------
    %Calculate the TCF with the newer numerical methods
    %----------------------------------------------------------------------
    if clockMode == 1
        tic;
    end
    % function C2 = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,timeArray)
    C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
    display_str = 'C2Maker 3state123 cyclical';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    C2_sim_plot = plot(C2_exp_x,C2_sim,'b--','LineWidth',2,'DisplayName',display_str);
    %RESULT:      Took 0.8411 seconds to run C2Maker 3state123 cyclical (192
    %times slower than the analytical method). For certain values this
    %fails/
    
    
    %Clean up the plot
    logx;
    xlabel('Time (sec)','FontSize',12);
    ylabel('C^{(2)}(\tau)','FontSize',12);
    legend('show');
    
    %Make a zeroline
    hold on
    line([C2_exp_x(1)+1e-6 C2_exp_x(end)],[0 0],'Color','red','linestyle',':','Linewidth',1);
    fprintf(['t12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
        '\n A1 = %f, A2 = %f, A3 = %f,\r\n'],...
        t12,t13,t21,t23,t31,t32,A1,A2,A3);
    
    
    %--------------------------------------------------------------------------
    % Save the data
    %--------------------------------------------------------------------------
    if saveMode == 1
        %--------------------------------------------------------------------------
        % Make a folder to hold output
        %-------------------------------------------------------------------------
        %Make output folder if it doesnt exist
        outputFolderName = [programName '_output'];
        if exist(outputFolderName,'dir') ~= 7
            mkdir(outputFolderName);
            disp('Making a folder to hold the output');
        end
        
        %--------------------------------------------------------------------------
        % Save a figure
        %-------------------------------------------------------------------------
        set(gcf,'PaperPositionMode','auto');
        extension = 'png';
        runNum = 1;
        foutName = ['C2simulationComparison_' num2str(runNum)];
        filePathOut = [outputFolderName filesep() foutName];
        while exist([filePathOut '.' extension],'file') == 2
            runNum = runNum + 1;
            fprintf('Increasing the run num to %f',runNum);
            foutName = ['C2simulationComparison_' num2str(runNum)];
            filePathOut = [outputFolderName filesep() foutName];
        end
        print(gcf,filePathOut,['-d' extension],'-r0');
        disp(['Saving the figure as ' foutName]);
    end
end
%//////////////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------
% Compare the 4-point TCFs
%--------------------------------------------------------------------------
%//////////////////////////////////////////////////////////////////////////
if compareC4Mode == 1
    figure(23);
    clf;
    set(gcf,'Color','w');
    set(gcf,'Name','C4: 3state123 Cyclical');
    
    disp('Comparing the time it takes to calculate 4-point TCFs...');
    
    %Make an array of time for the x-axis
    Npts = 150;
    time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
    
    tau1range = time;
    tau2 = 0;
    
    %----------------------------------------------------------------------
    %Calculate the C4 with the old analytical expressions (rewritten)
    %----------------------------------------------------------------------
    
    if clockMode == 1
        tic;
    end
    
    [C4_sim_ana,~,~] = C4maker_3state123_cyclical_analytical(t12,t13,t21,t23,t31,A1,A2,A3,tau2,tau1range);
    display_str = 'C4maker 3state123 cyclical analytical';
    
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    
    if plotMode == 1
        hold on;
        surf_TCF4point = surf(tau1range, tau1range', C4_sim_ana,'DisplayName',display_str);
        title('Four-point TCF: C^{(4)}','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('Time (\tau_3)','FontSize',14);
        zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
        
        view(28,36);
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
        
        drawnow();
        hold on;
    end
    
    %----------------------------------------------------------------------
    %Calculate the TCF with the newer numerical methods
    %----------------------------------------------------------------------
    if clockMode == 1
        tic;
    end
    
    % function [C4,C4_diff,C2] = C4Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,tau2,timeArray)
    [C4_sim_num,~,~] = C4Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,tau2,tau1range);
    display_str =  'C4Maker 3state123 cyclical (numerical)';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    %
    if plotMode == 1
        
        hold on;
        plot_TCF4point_analytical = mesh(tau1range, tau1range', C4_sim_num,'DisplayName',display_str);
        plot_TCF4point_analytical.FaceAlpha = 0.5;
        title('Four-point TCF: C^{(4)}','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('Time (\tau_3)','FontSize',14);
        zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
        
        view(28,36);
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
        lgd = legend;
        lgd.Location = "north";
        
    end
    
    
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    %--------------------------------------------------------------------------
    % Save the data
    %--------------------------------------------------------------------------
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if saveMode == 1
        %--------------------------------------------------------------------------
        % Make a folder to hold output
        %-------------------------------------------------------------------------
        %Make output folder if it doesnt exist
        outputFolderName = [programName '_output'];
        if exist(outputFolderName,'dir') ~= 7
            mkdir(outputFolderName);
            disp('Making a folder to hold the output');
        end
        
        %--------------------------------------------------------------------------
        % Save a figure
        %-------------------------------------------------------------------------
        set(gcf,'PaperPositionMode','auto');
        extension = 'png';
        runNum = 1;
        foutName = ['C4simulationComparison_' num2str(runNum)];
        filePathOut = [outputFolderName filesep() foutName];
        while exist([filePathOut '.' extension],'file') == 2
            runNum = runNum + 1;
            fprintf('Increasing the run num to %f',runNum);
            foutName = ['C4simulationComparison_' num2str(runNum)];
            filePathOut = [outputFolderName filesep() foutName];
        end
        print(gcf,filePathOut,['-d' extension],'-r0');
        disp(['Saving the figure as ' foutName]);
    end
    
end
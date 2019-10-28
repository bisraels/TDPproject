function statsComparison_3state123_linear(t12,t21,t23,t32,A1,A2,A3)
%--------------------------------------------------------------------------
%USER OPTIONS
%--------------------------------------------------------------------------
saveMode = 0;
clockMode = 1;%if clockMode == 1, tic; end if clockMode == 1, disp(['     Took ' num2str(toc) ' seconds to run ' display_str]); end
plotMode = 1;
verboseMode = 0 ;

compareHistMode = 1;
compareC2Mode = 1;
compareC4Mode = 1;

%--------------------------------------------------------------------------
% EXTERNAL PROGRAMS CALLS: (needs these codes in MATLAB PATH to function)
%--------------------------------------------------------------------------
% ANALYTICAL Algorithms (fast)
% (1) histMaker_3state123_linear_analytical    % Calculates Peq
% (2) C2maker_3state123_linear_analytical      % Calculates C2
% (3) C4maker_3state123_linear_analytical      % Calculates C4

% NUMERICAL Algorithms (slow)
% (0) ODEsolver_3state123_linear.m    % Calculates the conditional probabilities (symCondProb_3state123_linear.mat)
% (1) histMaker_3state123_linear      % Calculates Peq
% (2) C2Maker_3state123_linear        % Calculates C2
% (3) C4Maker_3state123_linear        % Calculates C4

programName = 'statsComparison_3state123_linear';
%--------------------------------------------------------------------------
% SET PARAMATERS
%--------------------------------------------------------------------------
%MODEL: 1 <--> 2 <--> 3
programName = 'C2Maker_3state123_linear.m';
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
        
end

%--------------------------------------------------------------------------
% Set the rates (4 rates)
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
% Compare the histograms
%--------------------------------------------------------------------------
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
    % Calculate the histogram using the analytical code
    %----------------------------------------------------------------------
    if clockMode == 1
        tic;
    end
    [Peq] = histMaker_3state123_linear_analytical(t12,t21,t23,t32,A1,A2,A3);
    display_str = 'histMaker_3state123_linear_analytical';
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
    fitHistPlot2.DisplayName = ['histMaker 3state123 linear analytical'];
    hold on;
    
    leg = legend;
    leg.Location = 'best';
    %----------------------------------------------------------------------
    % Calculate the histogram using the new code
    %----------------------------------------------------------------------
    hold on
    if clockMode == 1
        tic;
    end
    [Peq] = histMaker_3state123_linear(t12,t21,t23,t32,A1,A2,A3);
    display_str = 'histMaker_3state123_linear';
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
    fitHistPlot2.DisplayName = ['histMaker 3state123 linear'];
    
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



%--------------------------------------------------------------------------
% Compare the 2-point TCFs
%--------------------------------------------------------------------------
if compareC2Mode == 1
    disp('Comparing the time it takes to calculate 2-point TCFs...');
    figure(22)
    clf
    set(gcf,'Color','w');
    
    %Make an array of time for the x-axis
    Npts = 150;
    time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
    C2_exp_x = time;
    
    %     %----------------------------------------------------------------------
    %     %Calculate the TCF with the old analytical expressions (test)
    %     %----------------------------------------------------------------------
    %     if clockMode == 1
    %         tic;
    %     end
    %     %  tcf = C2maker_3state123_linear_analytical(t12,t13,t21,t23,t31,A1,A2,A3,time)
    %     tcf_simv2 = C2maker_3state123_linear_analytical_v2(t12,t21,t23,t32,A1,A2,A3,time);
    %     display_str = 'C2maker 3state123 linear analytical v2';
    %     if clockMode == 1
    %         disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    %     end
    %     tcf_plot = plot(C2_exp_x,tcf_simv2,'g:','LineWidth',2,'DisplayName',display_str);
    %     hold on;
    
    
    %----------------------------------------------------------------------
    %Calculate the TCF with the old analytical expressions (cleaner)
    %----------------------------------------------------------------------
    if clockMode == 1
        tic;
    end
    
    %  tcf = C2maker_3state123_linear_analytical(t12,t13,t21,t23,t31,A1,A2,A3,time)
    tcf_sim_ana = C2maker_3state123_linear_analytical(t12,t21,t23,t32,A1,A2,A3,time);
    display_str = 'C2maker 3state123 linear analytical';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    tcf_plot = plot(C2_exp_x,tcf_sim_ana,'r-','LineWidth',2,'DisplayName',display_str);
    hold on;
    %RESULT:      Took 0.0043612 seconds to run C2maker 3state123 linear analytical
    
    %----------------------------------------------------------------------
    %Calculate the TCF with the newer numerical methods
    %----------------------------------------------------------------------
    if clockMode == 1
        tic;
    end
    % function C2 = C2Maker_3state123_linear(t12,t13,t21,t23,t31,A1,A2,A3,timeArray)
    C2_sim = C2Maker_3state123_linear(t12,t21,t23,t32,A1,A2,A3,time);
    display_str = 'C2Maker 3state123 linear';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    C2_sim_plot = plot(C2_exp_x,C2_sim,'b--','LineWidth',2,'DisplayName',display_str);
    %RESULT:      Took 0.8411 seconds to run C2Maker 3state123 linear (192
    %times slower than the analytical method). For certain values this
    %fails/
    
    %----------------------------------------------------------------------
    %Calculate the TCF with the new numerical expressions (test)
    %----------------------------------------------------------------------
    % if clockMode == 1
    %     tic;
    % end
    % %  tcf = C2maker_3state123_linear_analytical(t12,t13,t21,t23,t31,A1,A2,A3,time)
    % tcf_numv2 = C2Maker_3state123_linear_v2(t12,t21,t23,t32,A1,A2,A3,time);
    % display_str = 'C2maker 3state123 linear v2';
    % if clockMode == 1
    %     disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    % end
    % tcf_plot = plot(C2_exp_x,tcf_numv2,'g:','LineWidth',4,'DisplayName',display_str);
    % hold on;
    
    
    %----------------------------------------------------------------------
    %Clean up the plot
    %----------------------------------------------------------------------
    logx;
    xlabel('Time (sec)','FontSize',12);
    ylabel('C^{(2)}(\tau)','FontSize',12);
    legend('show');
    
    %Make a zeroline
    hold on
    line([C2_exp_x(1)+1e-6 C2_exp_x(end)],[0 0],'Color','red','linestyle',':','Linewidth',1);
    fprintf(['t12 = %f, t21 = %f, t23 = %f, t32 = %f'...
        '\n A1 = %f, A2 = %f, A3 = %f,\r\n'],...
        t12,t21,t23,t32,A1,A2,A3);
    
    
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

%--------------------------------------------------------------------------
% Compare the 4-point TCFs
%--------------------------------------------------------------------------
if compareC4Mode == 1
    disp('Comparing the time it takes to calculate 4-point TCFs...');
    if plotMode == 1
        figure(23);
        clf;
        
        set(gcf,'Color','w');
        set(gcf,'Name','C4: 3state123 linear');
        
    end
    %Make an array of time for the x-axis
    Npts = 150;
    time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
    
    tau1range = time;
    tau2 = 0;
%     %----------------------------------------------------------------------
%     %Calculate the C4 with the old analytical expressions
%     %----------------------------------------------------------------------
%     if clockMode == 1
%         tic;
%     end
%     %     [kappa,theta,tcf] = FourPtTCF_3state(tau1range,tau2,tau3range,A0,A1,A2,k01,k10,k12,k21);
%     
%     disp('>>: Running FourPtTCF_3state_V2');
%     [~,~,C4_sim_ana1] = FourPtTCF_3state_V2(tau1range,tau2,tau1range',A1,A2,A3,k12,k21,k23,k32);
%     display_str = 'FourPtTCF_3state';
%     if clockMode == 1
%         disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
%     end
%     
%     if plotMode == 1
%         set(gcf,'Color','w');
%         set(gcf,'Name','C4');
%         
%         surf_C4_sim_ana1 = surf(tau1range, tau1range', C4_sim_ana1,'DisplayName','FourPtTCF 3state V2');
%         title(' Four-point TCF: C^{(4)}','FontSize',18)
%         xlabel('Time (\tau_1)','FontSize',14);
%         ylabel('Time (\tau_3)','FontSize',14);
%         zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
%         
%         view(28,36);
%         ax = gca;
%         ax.XScale = 'log';
%         ax.YScale = 'log';
%         
%         drawnow();
%         hold on;
%     end
    %----------------------------------------------------------------------
    %Calculate the C4 with the old analytical expressions (rewritten)
    %----------------------------------------------------------------------
    
    if clockMode == 1
        tic;
    end
    disp('>>: Running C4maker_3state123_linear_analytical');
    %     [C4,C4_diff,C2] = C4maker_3state123_linear_analytical(t12,t21,t23,t32,A1,A2,A3,tau2,tau1range)
    [C4_sim_ana,~,~] = C4maker_3state123_linear_analytical(t12,t21,t23,t32,A1,A2,A3,tau2,tau1range);
    display_str = 'C4maker 3state123 linear analytical';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    if plotMode == 1
        figure(23)
        C4_sim_ana = surf(tau1range, tau1range, C4_sim_ana,'DisplayName',display_str);
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
    disp('>>: Running C4Maker_3state123_linear');
    [C4_sim_Num,~,~] = C4Maker_3state123_linear(t12,t21,t23,t32,A1,A2,A3,tau2,tau1range);
    display_str =  'C4Maker 3state123 linear (numerical)';
    if clockMode == 1
        disp(['     Took ' num2str(toc) ' seconds to run ' display_str]);
    end
    
    if plotMode == 1
        figure(23)
        hold on;
        plot_C4_sim_Num = mesh(tau1range, tau1range, C4_sim_Num,'DisplayName',display_str);
        title('Four-point TCF: C^{(4)}','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('Time (\tau_3)','FontSize',14);
        zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
        
        view(28,36);
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
        legend('show')
        
            drawnow();
            hold on;
    end
    
    
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
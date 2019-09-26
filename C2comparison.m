% function C2comparison(t12,t13,t21,t23,t31,A1,A2,A3)
programName = 'C2comparison';
%--------------------------------------------------------------------------
%USER OPTIONS
%--------------------------------------------------------------------------
saveMode = 0;
compareHistMode = 1;
compareC2Mode = 1;
%--------------------------------------------------------------------------
% SET PARAMATERS
%--------------------------------------------------------------------------
% switch nargin
%     case 0
disp('Using default values in C2comparison');

% 
% t12_bounds = [1e-6,1000e-6];  %Paramater #1 is high--> med
% t13_bounds = [100e-6,10e-3];    %Paramater #2 is high --> low
% t21_bounds = [1e-6,1e-3];%Paramater #3 is med --> high
% t23_bounds = [1e-6,10e-3];%Paramater #4 is med --> low
% t31_bounds = [10e-6,10e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
% % *t32 wll be determined by the other rates
% 
% A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
% A2_bounds = [0.45,0.65];%Paramater #7 % Med FRET state
% A3_bounds = [0.30,0.45];%Paramater #8 %Low FRET state
% boundsArray = [t12_bounds;t13_bounds;t21_bounds;t23_bounds;t31_bounds;A1_bounds;A2_bounds;A3_bounds];
% 
% population = rand(1,8);
% for param_idx = 1:8
%     %To pick a random number in the interval of LB to UB:
%     % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
%     population(param_idx) = boundsArray(param_idx) + population(param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
% end
% t12 = population(1);
% t13 = population(2);
% t21 = population(3);
% t23 = population(4);
% t31 = population(5);
% A1 = population(6);
% A2 = population(7);
% A3 = population(8);

% end

k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;

k32 = k12*k23*k31/(k13*k21); % Detailed balance condition: %k31 will be the rate fixed by the others
t32 = 1/k32;

%--------------------------------------------------------------------------
% Compare the histograms
%--------------------------------------------------------------------------
if compareHistMode == 1
    figure(21);
    set(gcf,'Color','w');
  
    sigma_A1 = 0.15;
    sigma_A2 = 0.1;
    sigma_A3 = 0.1;
    FRET_bins = linspace(0,1,100);
    [p1_eq,p2_eq,p3_eq] = cyclic3state_hist(A1,A2,A3,k12,k21,k23,k32,k31);
    hist_sim = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim = sum(hist_sim);
    hist_sim = hist_sim./sum(hist_sim);
    fitHistPlot1 = plot(FRET_bins,hist_sim);
    fitHistPlot1.LineStyle = '-';
    fitHistPlot1.Color = 'red';
    fitHistPlot1.LineWidth = 2;
    fitHistPlot1.DisplayName = ['FourPtTCF_cyclic3state_xo'];
    
    hold on
    [Peq] = histMaker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3);
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
    fitHistPlot2.DisplayName = ['histMaker_3state123_cyclical'];
end
%--------------------------------------------------------------------------
% Compare the 2-point TCFs
%--------------------------------------------------------------------------
%Get ready for plotting
if compareC2Mode == 1
    figure(22)
    clf
    set(gcf,'Color','w');
    
    %Make an array of time for the x-axis
    Npts = 150;
    timeArray = [0:9,logspace(1,6.4771212,Npts)]/1e6;
    C2_exp_x = timeArray;
    
    %Calculate the TCF with the old analytical expressions
    % tcf = TCF_cyclic3state(time,A0,A1,A2,k01,k10,k12,k21,k20)
    tcf = TCF_cyclic3state(C2_exp_x,A1,A2,A3,k12,k21,k23,k32,k31);
    display_str = 'TCF cyclic3state';
    tcf_plot = plot(C2_exp_x,tcf,'r-','LineWidth',2,'DisplayName',display_str);
    hold on;
    
    %Calculate the TCF with the newer numerical methods
    C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
    display_str = 'C2Maker 3state123 cyclical';
    C2_sim_plot = plot(C2_exp_x,C2_sim,'b--','LineWidth',2,'DisplayName',display_str);
    
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
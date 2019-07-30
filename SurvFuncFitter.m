% function SurvFuncFitter(xIN,yIN,res)
%__________________________________________________________________________
% AUTHOR: Brett Israels & Claire Albrecht
%
% NAME: SurvFuncFitter.m
%
% FUNCTION: % Fits the Survival probability function created from dwell
% time distributions to a single exponential decay
%
% INPUT: 1) x: An array of times
%        2) y: The survival probability (1 --> 0)
%
% OUTPUT: 1) The time constant that charachterizes the rate
%         2) A figure which shows the data and the best fit
%         3) 
%
% MODIFICATION LOG:
clf
fileNameKeyword = 'survProb_';
fileNames = dir([fileNameKeyword '*']);
for fileName_idx = 1:numel(fileNames)
    clf;    %Clear the current Figure
    load(fileNames(fileName_idx).name,'time','prob','rate_ID')
 

functionName = 'SurvFuncFitter';

%----------------------------------------------------------------------
%  Set  paramaters
%----------------------------------------------------------------------

res = 0.03;
saveMode = 1;
saveFigMode = 1;
verboseMode = 0;
showProgressMode = 0;
findBestFit_mode = 1;
labelPlotWithFit_mode = 1;
%----------------------------------------------------------------------
%  Prepare for output
%----------------------------------------------------------------------
if saveMode
    outputFolderName = [functionName '_output' ];
    if exist(outputFolderName,'dir') ~= 7
        fprintf('     ***Making a folder called %s to store output.\r\n',outputFolderName);
        mkdir(outputFolderName);
    end
end

%----------------------------------------------------------------------
%  Fit the Survival Function to a single exponential
%----------------------------------------------------------------------
% Set up fittype and options.
[xData, yData] = prepareCurveData( time, prob );
fitFunction = 'exp(-x/t1)';

fprintf('Fitting a function to the raw data: y(x) = %s\r\n',fitFunction);

ft = fittype( fitFunction, 'independent', 'x', 'dependent', 'y' );

% Choose Fitting options
chosen_method = 'NonlinearLeastSquares';
opts = fitoptions( 'Method', chosen_method );
opts.Display = 'off';
fprintf('\t Fit options: %s\n',chosen_method);

%--------------------------------------------------------------------------
% Assign the Fit Paramaters: fit options
% LowerBound (LB), Upper Bound (UB) and StartingPoint (SP)
% fitFunction = 'exp(-x/t1)';
%--------------------------------------------------------------------------
t1_LB = res;
t1_UB = max(xData);
t1_SP = (t1_LB + t1_UB)/2;

%*** THESE MUT BE IN ALPHABETICAL ORDER ***
opts.Lower =      [t1_LB];
opts.Upper =      [t1_UB];

% Fit model to data.
randMultiplier = 1;
maxTrials = 10;
for i = 1:maxTrials
    t1_SPm = t1_SP + t1_SP*randn*randMultiplier;
    
    opts.StartPoint = [t1_SPm];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the Fit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [fitresult, gof] = fit( xData, yData, ft, opts );
    rsquare = gof.rsquare;
    rsquareNew = rsquare;
    if verboseMode
        disp(['i = ' num2str(i) '/' num2str(maxTrials)]);%
        disp(['    rsquare = ' num2str(rsquare)]);
        disp(fitresult);
    end
    
    if showProgressMode
        figure(10)
        plot(i,rsquare,'--gs',...
            'LineWidth',2,...
            'MarkerSize',10,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[0.5,0.5,0.5])
        ylim([.9,1]);
        hold on;
    end
    %-------------------------------------------------------------------------
    % Save the fit values
    %-------------------------------------------------------------------------
    
    t1 = fitresult.t1;
    t = 0:res:max(time);
    y = exp(-t/t1);
    
    if saveMode
        fitFunctionName = '1exp_fitresult';
        foutName_prefix = [rate_ID '_' fitFunctionName];
        filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
        if findBestFit_mode
            if exist(filepath_fit,'file') ~= 2
                disp('Did not find the file. Saving the initial fit');
                iterNumber = 1;
                save(filepath_fit,'fitresult','gof','xData','yData','t1','t','y','rsquare','iterNumber');
            else
                load(filepath_fit,'rsquare','iterNumber');
                iterNumber = iterNumber + 1 ;
                %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
                if rsquareNew > rsquare
                    disp(['*** Found a better fit: ' num2str(rsquareNew) ' > ' num2str(rsquare) '. Saving.']);
                    rsquare = rsquareNew;
                    save(filepath_fit,'fitresult','gof','xData','yData','t1','t','y','rsquare','iterNumber');
                else
                    %                         disp(['     appending the iteration number to' filepath]);
                    save(filepath_fit,'iterNumber','-append')
                end
            end
        else
            disp(['     Saving the fit results as ' foutName_prefix '.mat in ' outputFolderName]);
            save(filepath_fit,'fitresult','gof','xData','yData','t1','t','y','rsquare','iterNumber');
        end
    end
    fprintf('%d/%d fitresult #%d:  t1 = %f\n',i,maxTrials,iterNumber,t1);
end
%     save(filepath_fit,'xData','yData','opts','-append');
%-------------------------------------------------------------------------
% Plot the [Best] Fit with Data
%-------------------------------------------------------------------------
load(filepath_fit,'fitresult','gof','xData','yData','t1','t','y','rsquare','iterNumber');

t1 = fitresult.t1;
    t = 0:res:max(xData);
    y = exp(-t/t1);
    
fig = figure(1);
fig.Name = '1exp Fit';

plot(xData,yData,'b.','MarkerSize',10,'DisplayName','C^{(2)}(\tau)')
hold on;

load(filepath_fit,'fitresult','gof','t1','y','rsquare','iterNumber');
plot(t,y,'color','red','LineWidth',2,'DisplayName',['Best Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
hold on;

legend('show');
lgd = legend;
lgd.Location = 'NorthEast';
lgd.FontSize = 12;

% [sample_description,save_prefix] = sample_descriptionGetter()
sample_description = '3^{\prime}p(dT)_{15}';
title_str = [sample_description '{\color{blue}  rate ' rate_ID '}{ \color{red}Res = ' num2str(res) 'sec}'...
    10 'y = ' fitFunction];
title(title_str,'FontSize',16);
xlabel('time (sec)','fontsize',16);
ylabel('C^{(2)}(\tau)','fontsize',16);


% axis tight;
set(gcf,'Color',[1 1 1]);
grid on;
set(gca,'yscale','lin');
set(gca,'xscale','lin');
set(gca,'FontSize',14);
axis tight;
axis square

if labelPlotWithFit_mode
    
    text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1,'%2.3f') ' sec'],'Color','red','FontSize',14);
end

drawnow();

%-------------------------------------------------------------------------
% Save the Average TCF overlayed with the Exp Fit
%-------------------------------------------------------------------------
if saveFigMode
    fitFunctionName = '1exp_fitresult';
    
    foutName_prefix = [rate_ID '_' fitFunctionName];
        
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    extension = 'png';
    disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
    print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
    saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
end
end

%     plot_error_mode = 0;
%     if plot_error_mode
%         rectangle('Position',[x_begin,y_begin,x_width,y_height],'FaceColor',[0 .5 .5])
%     end
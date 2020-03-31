%--------------------------------------------------------------------------
% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  November 2019
%
% PURPOSE:  Create histograms for any model
%
% INPUT:    (1) Conditional Prpbability matrix (P, from k2P)
%           (2) The FRET values for the model  (A, a vector)
%
% OUTPUT:   (1) FRET Histogram
%
% MODIFICATION:
%
%--------------------------------------------------------------------------

% global plotMode verboseMode


% function [Peq, histPlot] = histMaker_Nstate(p, A, sigma_A)
function [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate(P, p, A, sigma_A)
global genNum plotMode verboseMode
%--------------------------------------------------------------------------
% PARAMETERS TO GET THE FUNCTION WORKING:
% switch nargin
% %     case 0
% plotMode = 1;
% verboseMode = 1;
% FRET_bins = 1:100;
% % A = sym('A',[N,1]);
% A = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7]; % example to get working
% syms t
% 
% K = [-25,  21,     31,         0,     0,       0,       0;...   %1
%     12,    -93,    2852/91,    42,    52,      0,       0;...   %2
%     13,    23,  -(31+2852/91), 0,     0,       0,       0;...   %3
%     0,     24,  0,             -42,   0,       0,       0;...   %4
%     0,     25,  0,              0,    -108,    65,      0;...   %5
%     0,     0,   0,              0,      56,  -(65+67),  76;...  %6
%     0,     0,   0,              0,       0,     67,     -(76)]; %7
% p = k2P(K,t);
% sigma_A1 = 0.15;
% sigma_A2 = 0.1;
% sigma_A3 = 0.1;
% sigma_A4 = 0.1;
% sigma_A5 = 0.1;
% sigma_A6 = 0.1;
% sigma_A7 = 0.1;
% sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4; sigma_A5; sigma_A6; sigma_A7];
% 
% clf;
% genNum = 1;
% end
% %--------------------------------------------------------------------------
%%

long_time = 3000000;
time = long_time;
t = time;
peq = diag(p);
Peq = double(subs(peq));

if verboseMode == 1
    if (sum(Peq) - 1)  < 1e-9
        disp('Equilibrium probabilities sum to 1! (within 1e-9 error)')
    else
        disp('Equilibrium probabilities do NOT sum to 1! (within 1e-9 error)')
        disp('The difference from one is:')
        disp(sum(Peq)  - 1)
        
    end
end

N = length(p);  % Number of states


% Need a list of Peq names for the legend of these histograms!
    Peq_sym = sym('Peq_',[N,1]);            % Define list of strings for the legend of the Peq fits
    Peq_char = [];
    Peq_name = [];
    Peq_names = [];
    for i = 1:N
        Peq_char = char(Peq_sym(i));        % Turn syms to char in a list
        Peq_name = [Peq_name; Peq_char];    % concatenate the list of chars
        Peq_names = [Peq_names; strrep(Peq_name(i,:),'eq','^{eq}')];  % Make the eq a super script
    end

%% SOMETHING IS WRONG WITH THE HIST_SIM_MAT CALCULATION!


clear hist_sim denom_hist_sim hist_sim_mat
% hist_sim_mat = [];
hist_sim = 0;
% hist_sim_mat = 0;
for i = 1:N
    hist_sim_temp = Peq(i) * exp(-(((FRET_bins - A(i))/sigma_A(i)).^2));
%     hist_sim_mat = [hist_sim_mat; hist_sim_temp];
    hist_sim = hist_sim + hist_sim_temp;
end

denom_hist_sim= sum(hist_sim(:));

disp('denom_hist_sim =')
disp(denom_hist_sim)        % This should be greater than 1.

% hist_sim = sum(hist_sim_mat)./denom_hist_sim;
hist_sim = hist_sim./denom_hist_sim;

disp('hist_sim shape =')
disp(size(hist_sim))

disp('max(hist_sim)=')
disp(max(hist_sim))

%%

if plotMode == 1
    figure(1);
    
    set(gcf,'Color','w');
    set(gcf,'Name','FRET Histogram');
    lineStyle = char('r--','g--','b--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
    
    FRET_bins = linspace(0,1,100);
    
    
    
    % Plot histogram of each state
    for i = 1:N
        %         plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2),...
        %             lineStyle(i,:),'LineWidth',1,'DisplayName',Peq_names(i,:));
        
        plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
            lineStyle(i,:),'LineWidth',1,'DisplayName',Peq_names(i,:));
        %         plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2)./sum(hist_sim),lineStyle(i,:),'LineWidth',1,'DisplayName',[Peq_names(i,:) '=' num2str(Peq(i))]);
        
        lgd = legend('show');
        %             lgd.Location = 'northwest';
        lgd.FontSize = 14;
        hold on
    end
    
    set(gca, 'FontSize', 14);
    
    
    if exist('fitHistPlot','var') == 1
        delete(fitHistPlot)
    end
    
    hold on

    % Plot overall histogram 
        histPlot = plot(FRET_bins, hist_sim,...
            lineStyle(i),'LineWidth',1,'DisplayName','Total Fit');
        legend('show')
        %     saveas(histPlot_state,sprintf('histPlot_state%d.mat',i))
    
    hold off
    
 
end
%%
% state1_HistPlot = plot(FRET_bins,p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2)./denom_hist_sim,...
%     'c--','LineWidth',1,'DisplayName',['p1_{eq} =' num2str(p1_eq)]);
% state2_HistPlot = plot(FRET_bins,p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2)./denom_hist_sim,...
%     'm--','LineWidth',1,'DisplayName',['p1_{eq} =' num2str(p2_eq)]);
% state3_HistPlot = plot(FRET_bins,p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2)./denom_hist_sim,...
%     'g--','LineWidth',1,'DisplayName',['p1_{eq} =' num2str(p3_eq)]);
% legend('show');


% end


% hist_sim = sum(hist_sim_mat);
% disp('sum(hist_sim) =')
% disp(sum(hist_sim))

% denom_hist_sim =    19.0355;

% hist_sim = hist_sim./denom_hist_sim;   % Normalize




    % Final fit curve
%     plot(FRET_bins,hist_sim,'r-','LineWidth',3,'DisplayName','Final Fit');
%     xlabel('FRET Efficiency','FontSize',14);
%     ylabel('Frequency','FontSize',14);
%     title('Simulated Histograms','FontSize',14);
%     hold on
    
    % Add histogram fit from Gen1 to fig1 - howis this different than the
    % one above? The final fit?
%     fitHistPlot = plot(FRET_bins,hist_sim);
%     fitHistPlot.LineStyle = '-';
%     fitHistPlot.Color = 'red';
%     fitHistPlot.LineWidth = 2;
%     fitHistPlot.DisplayName = ['Gen' num2str(genNum)];
%  



%     histCum=0 * FRET_bins;
%     for i = 1:N
%         histState_temp = Peq(i) * exp(-((FRET_bins-A(i))/sigma_A(i)).^2);
%         histCum = histCum + histState_temp;
%     end





    
    %
    %
    
    % if exist('state1_HistPlot','var') == 1
    %     delete(state1_HistPlot);
    %     delete(state2_HistPlot);
    %     delete(state3_HistPlot);
    % end
    
    
    % histPlot_name_sym = sym('histPlot_state',[N,1]);
    % histPlot_char = [];
    % histPlot_name = [];
    % for i = 1:N
    %     histPlot_char_temp = char(histPlot_name_sym(i));
    %     histPlot_char = [histPlot_char; histPlot_char_temp];
    % end
    
    % histPlot_state = sym('histPlot_state',[N,1]);
    
   
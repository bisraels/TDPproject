%--------------------------------------------------------------------------
% AUTHOR:   Claire Albrecht & Brett Israels
%
% CREATED:  November 2019
%
% PURPOSE:  Create histograms for any model
%
% INPUT:    (1) Conditional Prpbability matrix (P, from k2P)
%           (2) The FRET value s for the model  (A)
%
% OUTPUT:   (1) FRET Histogram
%
% MODIFICATION:
%
%--------------------------------------------------------------------------

% global plotMode verboseMode


function [Peq] = histMaker_Nstate(p, A, sigma_A)
%--------------------------------------------------------------------------
% PARAMETERS TO GET THE FUNCTION WORKING:

plotMode = 1;
verboseMode = 1;
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

%--------------------------------------------------------------------------


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

N = length(Peq);


if plotMode == 1
    figure(1);
    
    set(gcf,'Color','w');
    set(gcf,'Name','FRET Histogram');
    
    
    FRET_bins = linspace(0,1,100);
    
    hist_sim = 0;
    for i = 1:numel(Peq)
        hist_sim_temp = Peq(i) * exp(-((FRET_bins - A(i).^2)));
        hist_sim = hist_sim + hist_sim_temp;
    end
    denom_hist_sim= sum(hist_sim);
    
    hist_sim = hist_sim./denom_hist_sim;   % Normalize
    
    plot(FRET_bins,hist_sim,'r-','LineWidth',2,'DisplayName','Final Fit');
    xlabel('FRET Efficiency','FontSize',14);
    ylabel('Frequency','FontSize',14);
    title('Simulated Histograms','FontSize',14);
    

    Peq_sym = sym('Peq_',[N,1]);            % Define list of strings for the legend of the Peq fits
    Peq_char = [];
    Peq_name = [];
    Peq_names = [];
    for i = 1:N
        Peq_char = char(Peq_sym(i));        % Turn syms to char in a list
        Peq_name = [Peq_name; Peq_char];    % concatenate the list of chars
        Peq_names = [Peq_names; strrep(Peq_name(i,:),'eq','^{eq}')];  % Make the eq a super script
    end
    
    
    lineStyle = char('r--','g--','b--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
    
    for i = 1:N
        plot(FRET_bins, Peq(i) * exp(-((FRET_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,lineStyle(i,:),'LineWidth',1,'DisplayName',Peq_names(i,:));
        lgd = legend('show');
        lgd.Location = 'northwest';
        lgd.FontSize = 14;
        hold on
    end
    hold off
    set(gca, 'FontSize', 14);
    
end

end


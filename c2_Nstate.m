% function [C2] = c2_nstate(A, K, tau1, tau2, tau3)

verboseMode = 1;

K = [-25,  21,     31,         0,     0,       0,       0;...   %1
    12,    -93,    2852/91,    42,    52,      0,       0;...   %2
    13,    23,  -(31+2852/91), 0,     0,       0,       0;...   %3
    0,     24,  0,             -42,   0,       0,       0;...   %4
    0,     25,  0,              0,    -108,    65,      0;...   %5
    0,     0,   0,              0,      56,  -(65+67),  76;...  %6
    0,     0,   0,              0,       0,     67,     -(76)]; %7
A = [.1:.1:.7];
N = length(K);



P_tau1 = k2P(K,tau1);
P_tau2 = k2P(K,tau2);
P_tau3 = k2P(K,tau3);

long_time = 300;
Peq = diag(k2P(K, long_time));

if verboseMode == 1
    if (sum(Peq) - 1) < 1e-9
        disp('Equilibrium probabilities sum to 1!')
    end
end

%--------------------------------------------------------------------------
% Calculate C2 using the 3x3xnumel(time) sized conditional probability
% Matrix P where P(j,i,:) is the prob of going from i->j in time t
%--------------------------------------------------------------------------
tic
Npts = numel(time);
C2_sim = zeros(size(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        C2_sim_temp = A(j) * reshape(P(j,i,:),[1,Npts]) * A(i) * Peq(i);
        C2_sim = C2_sim + C2_sim_temp;
    end
end
elapsedTime = toc;
disp(['Took ' num2str(elapsedTime) ' seconds to calculate C2 for N = ' num2str(N)]);
%Takes about 2 msec

%%
%--------------------------------------------------------------------------
%  Plot two point TCF
%--------------------------------------------------------------------------
% if plotMode == 1
    figure(2)
    clf;
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    TCF2pt = plot(time,C2_sim,'LineWidth',2,'DisplayName','C2');
    
    title('Two point TCF','FontSize',18)
    xlabel('Time','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);
    
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
    
% end


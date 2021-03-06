function [time, C2, C4] = P2C(P, K, time, A)
% FUNCTION: Computes the Correlation functions from the Conditional
% Probabilities
%MODIFICATIONS
% 20191216 BI Adding the FRET values as an input
% 20191217 BI Making sure Amean is calculated to be a number not a vector
%--------------------------------------------------------------------------

% global plotMode
global fitC2Mode fitC4Mode
clockMode = 0;
% User Options
plotMode = 0;

switch nargin
    case 0
        %Give a K matrix
        %     8state
        K = [...
            -(12 + 13), 21, 31, 0, 0, 0, 0, 0;...
            12, -(21 + 23 + 24 + 25), (12*23*31/(13*21)), 42, 52, 0, 0, 0;...
            13, 23, -(31 + ((12*23*31/(13*21)))), 0, 0, 0, 0, 0;...
            0, 24, 0, -42, 0, 0, 0, 0;...
            0, 25, 0, 0, -(52 + 56), 65, 0, 0;...
            0, 0, 0, 0, 56, -(65 + 67 + 68), 76, 86;...
            0, 0, 0, 0, 0, 67, -(76+78), 87;...
            0, 0, 0, 0, 0, 68, 78, -(86 + 87);...
            ];
        %Create the matrix of conditional probabilities
        %         [P, V, p, time] = k2P();
        Npts = 150;
        time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;
        
        [P, ~, ~, time] = k2P(K,time);
        
        [N,~,~] = size(P);
        A = linspace(.1,.9,N);
        
        clockMode = 1;
        fitC2Mode = 1;
        fitC4Mode = 1;
        plotMode = 1;
end

%Make a time array
% Npts = 150;
% time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;


% P = k2P(K,time);

%[N,~,timesteps] = size(P);
[N,~,~] = size(P);
timesteps = length(time);
% Make a row vector of the discrete probability
% x = diag(A) returns a column vector of the main diagonal elements of A.
Peq = diag(P(:,:,end))';

%Make a  row vector of the FRET Values
% A = linspace(.1,.9,N);

%Subtract off the mean of A
Amean = sum(A.*Peq);
A = A - Amean;

%--------------------------------------------------------------------------
% Calculate 2 point TCF (Sums using loops)
%--------------------------------------------------------------------------
if fitC2Mode == 1
    if exist('C2','var')
        clear('C2');
    end
    % syms C2;
    % C2 = [];
    if clockMode == 1
        tic
    end
    C2 = zeros(1,timesteps);
    for i = 1:numel(A)
        for j = 1:numel(A)
            C2temp = A(j)*reshape(P(j,i,:),[1 timesteps])*A(i)*Peq(i);
            C2 =  C2 + C2temp;
            
        end
    end
    if clockMode == 1
        elapsedTime =toc;
        disp(['Takes '  num2str(elapsedTime) ' seconds to calculate C2 for an ' num2str(N) ' state model']);
    end
else
    C2 = 0;
end
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


%--------------------------------------------------------------------------
% Calculate 2 point TCF (Matrix)
%--------------------------------------------------------------------------
%
% Amat = A.*eye(numel(A));
% % C2test = Amat*p*Amat*Peq;
% C2test = sum(Amat*cP*Amat*Peq);

% return
%--------------------------------------------------------------------------
% Calculate 4 point TCF
%--------------------------------------------------------------------------
if fitC4Mode == 1
    
    [P_tau2eq0] = k2P(K,0);
    C4 = zeros(timesteps,timesteps);
    if clockMode == 1
        tic
    end
    for i = 1:numel(A)
        for j = 1:numel(A)
            for k = 1:numel(A)
                for l = 1:numel(A)
                    
                    C4temp = A(l)*reshape(P(l,k,:),[1 timesteps])'*A(k)*P_tau2eq0(k,j)*A(j)*reshape(P(j,i,:),[1 timesteps])*A(i)*Peq(i);
                    C4 =  C4 + C4temp;
                    
                end
            end
        end
    end
    if clockMode == 1
        elapsedTime =toc;
        disp(['Takes '  num2str(elapsedTime) ' seconds to calculate C4 for an ' num2str(N) ' state model']);
    end
else
    C4 = 0;
end
if plotMode == 1
    figure(3)
    
    if exist('histPlot','var')  == 1
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

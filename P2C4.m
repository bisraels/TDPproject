function [C4,time] = P2C4(K,P,A,time,tau2)
clockMode = 0;
plotMode = 0;
switch nargin
    case 0
        
        [K,A,time,~] = paramSim_8state();
        P = K2P(K);
        tau2 = 0;
        
end
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
Ams = A - Amean;

%--------------------------------------------------------------------------
% Calculate 4point TCF (Sums using loops)
%--------------------------------------------------------------------------
if clockMode == 1
    tic
end
[P_tau2eq0] = k2P(K,tau2);
if clockMode == 1
    elapsedTime = toc;
    disp(['Time to calculate calculate Pji(t=0) for an ' num2str(N) ' state model: ' num2str(elapsedTime) ]);
end
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
    elapsedTime = toc;
    disp(['Time to calculate calculate C4 for an ' num2str(N) ' state model: ' num2str(elapsedTime) ]);
end

if plotMode == 1
    
    
    figure(3);
    set(gcf,'Color','w');
    set(gcf,'Name','4-point Time Correlation Function');
    
    surf(time,time,C4);
    shading interp;
    title_str = ['4-point TCF: \tau_2 = ' num2str(tau2)];
    title(title_str,'FontSize',12);
    xlabel('\tau_1 (sec)','fontsize',12);
    zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    grid on
    axis tight;
end



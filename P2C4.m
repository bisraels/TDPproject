function [C4,time] = P2C4(K,P,A,time,tau2)
clockMode = 0;
plotMode = 0;
switch nargin
    case 0
        
        [K,A,time,~] = paramSim_8state();
        P = K2P(K);
        
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
% Calculate 2 point TCF (Sums using loops)
%--------------------------------------------------------------------------
if clockMode == 1
    tic
end
C4 = zeros(timesteps,timesteps);
for i = 1:numel(Ams)
    for j = 1:numel(Ams)
        C2temp = Ams(j)*reshape(P(j,i,:),[1 timesteps])*Ams(i)*Peq(i);
        C2 =  C2 + C2temp;
    end
end
if clockMode == 1
    elapsedTime =toc;
    disp(['Takes '  num2str(elapsedTime) ' seconds to calculate C2 for an ' num2str(N) ' state model']);
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
% Calculate 2 point TCF (with matrix)
%--------------------------------------------------------------------------
%%
% A = linspace(0,1,N);
% Amat = A.*eye(N);
% Amat3D = ones(N,N,length(time));
% Peq2D = ones(N,length(time));
% for i = 1:length(time)
%     Amat3D(:,:,i) = Amat;
%     Peq2D(:,i) = Peq;
% end
% % C2test = sum(Amat*P*Amat*Peq);
% C2test = sum(Amat3D*P*Amat3D*Peq2D);
%
%


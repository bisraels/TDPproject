syms t tau1 tau3


time =  time_sim;
[P, ~, p, time] = k2P(K,time);

NumStates  = length(K);

% Make a column vector of the discrete probability
Peq = diag(P(:,:,end));
%Make a  row vector of the FRET Values
A = linspace(.1,.9,NumStates);

%Subtract off the mean of A
Amean = A*Peq;
A = A - Amean;


t  = tau1;
p_tau1 = subs(p);

t  = tau3;
p_tau3 = subs(p);


dp_tau1 = diff(p_tau1, tau1) * tau1;

dp_tau3 = diff(p_tau3, tau3) *  tau3;


tau1 = time;
tau3 = tau1;
timeSteps = length(tau1);

dP_tau1 = double(subs(dp_tau1(:)));
dP_tau1 = reshape(dP_tau1, [NumStates,NumStates,timeSteps]);

dP_tau3 = double(subs(dp_tau3(:)));
dP_tau3 = reshape(dP_tau3, [NumStates,NumStates,timeSteps]);


tic
D2 = zeros(1,timeSteps);
for i = 1:numel(A)
    for j = 1:numel(A)
            D2temp = A(j)*reshape(dP_tau1(j,i,:),[1 timeSteps])*A(i)*Peq(i);
            D2 =  D2 + D2temp;

    end
end
elapsedTime =toc;
disp(['Takes '  num2str(elapsedTime) ' seconds to calculate D2 for an ' num2str(NumStates) ' state model']);

if plotMode == 1
    figure(2)
    set(gcf,'Color','w');
    
    plot(tau1,D2);
    title_str = ['Two point time decay correlation function'];
    title(title_str,'FontSize',18);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('D^{(2)}(\tau)','fontsize',16);
    set(gca,'yscale','linear');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    grid on
    axis tight;
    
    drawnow();
end

%%

[P_tau2eq0] = k2P(K,0);
D4 = zeros(timeSteps,timeSteps);
tic
for i = 1:numel(A)
    for j = 1:numel(A)
        for k = 1:numel(A)
            for l = 1:numel(A)
                
                D4temp = A(l)*reshape(dP_tau3(l,k,:),[1 timeSteps])'*A(k)*P_tau2eq0(k,j)*A(j)*reshape(dP_tau1(j,i,:),[1 timeSteps])*A(i)*Peq(i);
                D4 =  D4 + D4temp;
                
            end
        end
    end
end
% D4 = double(D4);
elapsedTime =toc;
disp(['Takes '  num2str(elapsedTime) ' seconds to calculate D4 for an ' num2str(NumStates) ' state model']);

if plotMode == 1
    figure(4)
    set(gcf,'Color','w');
    
    surf(tau1,tau3,D4);
    title_str = ['Four point time decay function'];
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

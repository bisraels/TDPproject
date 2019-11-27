syms t tau1 tau3

[~, ~, p, ~] = k2P(K,time);


t  = tau1;
p_tau1 = subs(p);

t  = tau3;
p_tau3 = subs(p);


dp_tau1 = diff(p_tau1, tau1) * tau1;

dp_tau3 = diff(p_tau3, tau3) *  tau3;


tau1 = 1:100;
tau3 = tau1;
timeSteps = length(tau1);

dP_tau1 = double(subs(dp_tau1(:)));
dP_tau1 = reshape(dP_tau1, [N,N,timeSteps]);

dP_tau3 = double(subs(dp_tau3(:)));
dP_tau3 = reshape(dP_tau3, [N,N,timeSteps]);


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
disp(['Takes '  num2str(elapsedTime) ' seconds to calculate D4 for an ' num2str(N) ' state model']);

if plotMode == 1
    figure(4)
    set(gcf,'Color','w');
    
    surf(tau1,tau3,D4);
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

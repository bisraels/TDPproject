function C2comparison(t12,t13,t21,t23,t31,A1,A2,A3)
switch nargin
    case 0
        disp('Using default values');
        t12 = 1e-4;
        t13 = 0.0061;
        t21 = 3.27e-5;
        t23 = 1e-6;
        t31 = 1e-3;
        A1 = 0.7786;
        A2 = 0.6161;
        A3 = 0.4811;
        
end

k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;

k32 = k12*k23*k31/(k13*k21); % Detailed balance condition: %k31 will be the rate fixed by the others

%Get ready for plotting
figure(21)
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
C2plot = plot(C2_exp_x,tcf,'r-','LineWidth',2,'DisplayName',display_str);
hold on;

%Calculate the TCF with the newer numerical methods
C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
display_str = 'C2Maker 3state123 cyclical';
C2plot = plot(C2_exp_x,C2_sim,'b--','LineWidth',2,'DisplayName',display_str);

%Clean up the plot
logx;
logy;
xlabel('Time (sec)','FontSize',12);
ylabel('C^{(2)}(\tau)','FontSize',12);
legend('show');
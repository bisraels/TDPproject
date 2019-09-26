function C2comparison(t12,t13,t21,t23,t31,A1,A2,A3)
switch nargin
    case 0
        disp('Using default values in C2comparison');
        
        
t12_bounds = [1e-6,1000e-6];  %Paramater #1 is high--> med
   t13_bounds = [100e-6,10e-3];    %Paramater #2 is high --> low
t21_bounds = [1e-6,1e-3];%Paramater #3 is med --> high
t23_bounds = [1e-6,10e-3];%Paramater #4 is med --> low
t31_bounds = [10e-6,10e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
% *t32 wll be determined by the other rates

A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
A2_bounds = [0.45,0.65];%Paramater #7 % Med FRET state
A3_bounds = [0.30,0.45];%Paramater #8 %Low FRET state
boundsArray = [t12_bounds;t13_bounds;t21_bounds;t23_bounds;t31_bounds;A1_bounds;A2_bounds;A3_bounds];

population = rand(1,8);
for param_idx = 1:8
    %To pick a random number in the interval of LB to UB:
    % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
    population(param_idx) = boundsArray(param_idx) + population(param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
end
         t12 = population(1);
    t13 = population(2);
    t21 = population(3);
    t23 = population(4);
    t31 = population(5);
    A1 = population(6);
    A2 = population(7);
    A3 = population(8);
      
end

k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;

k32 = k12*k23*k31/(k13*k21); % Detailed balance condition: %k31 will be the rate fixed by the others


%--------------------------------------------------------------------------
% Compare the histograms
%--------------------------------------------------------------------------
% figure(21)
% set(gcf,'Color','w');


%--------------------------------------------------------------------------
% Compare the 2-point TCFs
%--------------------------------------------------------------------------
%Get ready for plotting
figure(22)
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
tcf_plot = plot(C2_exp_x,tcf,'r-','LineWidth',2,'DisplayName',display_str);
hold on;

%Calculate the TCF with the newer numerical methods
C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,C2_exp_x);
display_str = 'C2Maker 3state123 cyclical';
C2_sim_plot = plot(C2_exp_x,C2_sim,'b--','LineWidth',2,'DisplayName',display_str);

%Clean up the plot
logx;
xlabel('Time (sec)','FontSize',12);
ylabel('C^{(2)}(\tau)','FontSize',12);
legend('show');

line([C2_exp_x(1) C2_exp_x(end)],[0 0],'Color','red','linestyle',':','Linewidth',1);

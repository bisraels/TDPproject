function C2_sim = C2maker_3state123_linear_analytical_v2(t12,t21,t23,t32,A1,A2,A3,time)
verboseMode = 0;
plotMode = 1;
saveMode = 0;
%MODEL: 1 <--> 2 <--> 3
programName = 'C2Maker_3state123_linear.m';
switch nargin
    case 0
        disp(['Using default values in ' programName]);
        
        t12_bounds = [1e-6,1000e-6];  %Paramater #1 is high--> med
        t21_bounds = [1e-6,1e-3];%Paramater #3 is med --> high
        t23_bounds = [1e-6,10e-3];%Paramater #4 is med --> low
        t32_bounds = [10e-6,10e-3];  %Paramater #5 is low --> Medium
        
        A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
        A2_bounds = [0.45,0.65];%Paramater #7 % Med FRET state
        A3_bounds = [0.30,0.45];%Paramater #8 %Low FRET state
        boundsArray = [t12_bounds;t21_bounds;t23_bounds;t32_bounds;A1_bounds;A2_bounds;A3_bounds];
        
        Nparams = length(boundsArray);
        population = rand(1,Nparams);
        for param_idx = 1:Nparams
            %To pick a random number in the interval of LB to UB:
            % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
            population(param_idx) = boundsArray(param_idx) + population(param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
        end
        t12 = population(1);
        t21 = population(2);
        t23 = population(3);
        t32 = population(4);
        A1 = population(5);
        A2 = population(6);
        A3 = population(7);
        
        
        Npts = 150;
        time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
        
    case 7
        Npts = 150;
        time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
        
end
%--------------------------------------------------------------------------
% Set the rates (4rates)
%--------------------------------------------------------------------------
k12 = 1/t12;
k21 = 1/t21;
k23 = 1/t23;
k32 = 1/t32;

% Define eigenvalues lam1 and lam2 (3state123 Linear Analytical)
c1 = 0.5*(k12+k21+k23+k32);
c2 = 0.5*sqrt(k12^2 + 2*k12*k21 - 2*k12*k23 - 2*k12*k32 + k21^2 + 2*k21*k23 - 2*k21*k32 + k23^2 + 2*k23*k32 + k32^2);
lam1 = c1 + c2;
lam2 = c1 - c2;

% % SOLUTION FROM ODEsolver_3state123_linear
% % Eigenvalue 1 = 0
% % Eigenvalue 2 = - k12/2 - k21/2 - k23/2 - k32/2 - (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2
% % Eigenvalue 3 = (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2 - k21/2 - k23/2 - k32/2 - k12/2
% If c1 = k12/2 + k21/2 + k23/3 + k32/2
% then Eigenvalue 2 = -c1 +


if verboseMode == 1
    disp(['The final eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/lam1) ' microseconds '...
        ' tau2 = 1/eval2 = ' num2str(1e6*1/lam2) ' microseconds.']);
end

% % SOLUTION FROM ODEsolver_3state123_linear
% % Eigenvector 1 corresponds to the 0 eigenvalue
% % Eigenvector 1 = (k21*k32)/(k12*k23)
% % Eigenvector 1 = k32/k23
% % Eigenvector 1 = 1

% % Eigenvector 2 corresponds to the 1st non-zero eigenvalue
% % Eigenvector 2 = (k12/2 + k21/2 + k23/2 + k32/2 + (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23 - (k23 + k32)/k23
% % Eigenvector 2 = k32/k23 - (k12/2 + k21/2 + k23/2 + k32/2 + (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23
% % Eigenvector 2 = 1

% % Eigenvector 3 corresponds to the 3rd eigenvalue eigenvalue
% % Eigenvector 3 = (k12/2 + k21/2 + k23/2 + k32/2 - (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23 - (k23 + k32)/k23
% % Eigenvector 3 = k32/k23 - (k12/2 + k21/2 + k23/2 + k32/2 - (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23
% % Eigenvector 3 = 1

% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #. (original)
% % % a = (c1 + c2)/k21 - (k12 + k21)/k21;
% % % b = k12/k21 - (c1 + c2)/k21;
% % % c = 1;
% % %
% % % % Define eigenvector components corresponding to lam2: v2 = [x,y,z]. (original)
% % % x = (c1 - c2)/k21 - (k12 + k21)/k21;
% % % y = k12/k21 - (c1 - c2)/k21;
% % % z = 1;
% % Eigenvector 2 = (k12/2 + k21/2 + k23/2 + k32/2 + (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23 - (k23 + k32)/k23
% % Eigenvector 2 = (c1+c2)/k23 - (k23 + k32)/k23
% % Eigenvector 2 = k32/k23 - (k12/2 + k21/2 + k23/2 + k32/2 + (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23
% % Eigenvector 2 = 1

% Define eigenvector components corresponding to lam2: v2 = [x,y,z]. (altered)
a = (c1 + c2)/k23 - (k23 + k32)/k23;
b = k32/k23 - (c1 + c2)/k23;
c = 1;

% % Eigenvector 3 = (k12/2 + k21/2 + k23/2 + k32/2 - (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23 - (k23 + k32)/k23
% % x = (c1 - c2)/k23 - (k32 + k23)/k23;
% % Eigenvector 3 = k32/k23 - (k12/2 + k21/2 + k23/2 + k32/2 - (2*k12*k21 - 2*k12*k23 - 2*k12*k32 + 2*k21*k23 - 2*k21*k32 + 2*k23*k32 + k12^2 + k21^2 + k23^2 + k32^2)^(1/2)/2)/k23
% % Eigenvector 3 = 1
% Define eigenvector components corresponding to lam2: v2 = [x,y,z]. (altered)
x = (c1 - c2)/k23 - (k23 + k32)/k23;
y = k32/k23 - (c1 - c2)/k23;
z = 1;

% % % Define equilibrium values of populations (original)
% % p2_eq = 1/(k21/k12 + 1 + k23/k32);
% % p1_eq = k21/k12*p2_eq;
% % p3_eq = k23/k32*p2_eq;

% Define equilibrium values of populations (altered)
p2_eq = 1/(k23/k32 + 1 + k21/k12);
p3_eq = k23/k32*p2_eq;
p1_eq = k21/k12*p2_eq;

% % % ???Original
% % % Define constants needed for differential equation solutions, for various
% % % initial values of the population. m_2 is the constant "m" if the
% % % population of state 2 begins at 1 (others would then be zero).
% % n_3 = (-p2_eq - b/a*(1 - p3_eq))/(y - b/a*x);
% % n_2 = (1 - p2_eq + b/a*p3_eq)/(y - b/a*x);
% % n_1 = (-p2_eq + b/a*p3_eq)/(y - b/a*x);
% % m_3 = (1 - p3_eq - n_3*x)/a;
% % m_2 = (-p3_eq - n_2*x)/a;
% % m_1 = (-p3_eq - n_1*x)/a;
% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_3 = (-p2_eq - b/a*(1 - p1_eq))/(y - b/a*x);
n_2 = (1 - p2_eq + b/a*p1_eq)/(y - b/a*x);
n_1 = (-p2_eq + b/a*p1_eq)/(y - b/a*x);
m_3 = (1 - p1_eq - n_3*x)/a;
m_2 = (-p1_eq - n_2*x)/a;
m_1 = (-p1_eq - n_1*x)/a;

% Define populations as a function of time for different initial values.
% p2_1 population 2 as a function of time, with the population initially in
% state 1.
p3_3 = p3_eq + m_3*a*exp(-lam1*time) + n_3*x*exp(-lam2*time);
p3_2 = p3_eq + m_2*a*exp(-lam1*time) + n_2*x*exp(-lam2*time);
p3_1 = p3_eq + m_1*a*exp(-lam1*time) + n_1*x*exp(-lam2*time);
p2_3 = p2_eq + m_3*b*exp(-lam1*time) + n_3*y*exp(-lam2*time);
p2_2 = p2_eq + m_2*b*exp(-lam1*time) + n_2*y*exp(-lam2*time);
p2_1 = p2_eq + m_1*b*exp(-lam1*time) + n_1*y*exp(-lam2*time);
p1_3 = p1_eq + m_3*c*exp(-lam1*time) + n_3*z*exp(-lam2*time);
p1_2 = p1_eq + m_2*c*exp(-lam1*time) + n_2*z*exp(-lam2*time);
p1_1 = p1_eq + m_1*c*exp(-lam1*time) + n_1*z*exp(-lam2*time);

% Subtract mean values
Amean = p1_eq*A1 + p2_eq*A2 + p3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

% Calculate TCF
% Calculate C2_sim
C2_sim = p3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
    p2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
    p1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);

if plotMode == 1
    figure(2)
    
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    TCF2pt = plot(time,C2_sim,'LineWidth',2);
    
    title('Analytical Two point TCF','FontSize',18)
    xlabel('Time (\tau_1)','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);
    %xlim([10^-5 10^0])
    % ylim([C2(500) C2(0)])
    
    ax = gca;
    ax.XScale = 'log';
    set(gca,'yscale','log')
    
    if saveMode == 1
        saveName = ['C2_','example'];
        saveas(TCF2pt,saveName, 'png')
    end
end

function [C4,C4_diff,C2] = C4maker_3state123_linear_analytical(t12,t21,t23,t32,A1,A2,A3,tau2,tau1range)
%MODEL: 1 <--> 2 <--> 3
programName = 'C4maker_3state123_linear_analytical.m';
% disp([':>> Running ' programName '.m']);
verboseMode = 0;
plotMode = 0;
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
        
        tau2 = 0;
        
        Npts = 150;
        tau1range = [0:9,logspace(1,6.4771212,Npts)]/1e6;
        
    case 7
        tau2 = 0;
        Npts = 150;
        tau1range = [0:9,logspace(1,6.4771212,Npts)]/1e6;
end

%--------------------------------------------------------------------------
% Set the rates (4rates)
%--------------------------------------------------------------------------
k12 = 1/t12;
k21 = 1/t21;
k23 = 1/t23;
k32 = 1/t32;

%--------------------------------------------------------------------------
% Ready the time arrays
%--------------------------------------------------------------------------
tau3range = tau1range';

% Initialize 4-pt. TCF (C4) matrices
C4 = zeros(length(tau1range),length(tau3range));
C4_diff = zeros(length(tau1range),length(tau3range));
C2Product = zeros(length(tau1range),length(tau3range));

% Define eigenvalues lam1 and lam2 (3state123 Linear Analytical)
c1 = 0.5*(k12+k21+k23+k32);
c2 = 0.5*sqrt(k12^2 + 2*k12*k21 - 2*k12*k23 - 2*k12*k32 + k21^2 + 2*k21*k23 - 2*k21*k32 + k23^2 + 2*k23*k32 + k32^2);
lam1 = c1 + c2;
lam2 = c1 - c2;

if verboseMode == 1
    disp(['   ***The eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/lam1) ' microseconds '...
        ' tau2 = 1/eval2 = ' num2str(1e6*1/lam2) ' microseconds.']);
end


% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
a = (c1 + c2)/k21 - (k12 + k21)/k21;
b = k12/k21 - (c1 + c2)/k21;
c = 1;

% Define eigenvector components corresponding to lam2: v2 = [x,y,z].
x = (c1 - c2)/k21 - (k12 + k21)/k21;
y = k12/k21 - (c1 - c2)/k21;
z = 1;

% Define equilibrium values of populations
p2_eq = 1/(k21/k12 + 1 + k23/k32);
p1_eq = k21/k12*p2_eq;
p3_eq = k23/k32*p2_eq;

% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_3 = (-p2_eq - b/a*(1 - p3_eq))/(y - b/a*x);
n_2 = (1 - p2_eq + b/a*p3_eq)/(y - b/a*x);
n_1 = (-p2_eq + b/a*p3_eq)/(y - b/a*x);
m_3 = (1 - p3_eq - n_3*x)/a;
m_2 = (-p3_eq - n_2*x)/a;
m_1 = (-p3_eq - n_1*x)/a;



% % % Define populations as a function of time for different initial values.
% % % p2_1 population 2 as a function of time, with the population initially in
% % % state 1.
% % p3_3 = p3_eq + m_3*a*exp(-lam1*time) + n_3*x*exp(-lam2*time);
% % p3_2 = p3_eq + m_2*a*exp(-lam1*time) + n_2*x*exp(-lam2*time);
% % p3_1 = p3_eq + m_1*a*exp(-lam1*time) + n_1*x*exp(-lam2*time);
% % p2_3 = p2_eq + m_3*b*exp(-lam1*time) + n_3*y*exp(-lam2*time);
% % p2_2 = p2_eq + m_2*b*exp(-lam1*time) + n_2*y*exp(-lam2*time);
% % p2_1 = p2_eq + m_1*b*exp(-lam1*time) + n_1*y*exp(-lam2*time);
% % p1_3 = p1_eq + m_3*c*exp(-lam1*time) + n_3*z*exp(-lam2*time);
% % p1_2 = p1_eq + m_2*c*exp(-lam1*time) + n_2*z*exp(-lam2*time);
% % p1_1 = p1_eq + m_1*c*exp(-lam1*time) + n_1*z*exp(-lam2*time);

% Subtract mean values
Amean = p1_eq*A1 + p2_eq*A2 + p3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

% % % Calculate TCF
% % % Calculate C2_sim
% % C2_sim = p3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
% %          p2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
% %          p1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
% %      
t1 = tau1range;
t2 = tau2;
t3 = tau3range;
% disp(['Size t1 = ' num2str(size(t1))]);
% disp(['Size t3 = ' num2str(size(t3))]);

% Calculate TCF
path1 = p3_eq*A3*p3(3,t1)*A3*p3(3,t2)*A3.*p3(3,t3)*A3;
path2 = p3_eq*A3*p3(3,t1)*A3*p3(3,t2)*A3.*p2(3,t3)*A2;
path3 = p3_eq*A3*p3(3,t1)*A3*p3(3,t2)*A3.*p1(3,t3)*A1;
path4 = p3_eq*A3*p3(3,t1)*A3*p2(3,t2)*A2.*p3(2,t3)*A3;
path5 = p3_eq*A3*p3(3,t1)*A3*p2(3,t2)*A2.*p2(2,t3)*A2;
path6 = p3_eq*A3*p3(3,t1)*A3*p2(3,t2)*A2.*p1(2,t3)*A1;
path7 = p3_eq*A3*p3(3,t1)*A3*p1(3,t2)*A1.*p3(1,t3)*A3;
path8 = p3_eq*A3*p3(3,t1)*A3*p1(3,t2)*A1.*p2(1,t3)*A2;
path9 = p3_eq*A3*p3(3,t1)*A3*p1(3,t2)*A1.*p1(1,t3)*A1;

path10 = p3_eq*A3*p2(3,t1)*A2*p3(2,t2)*A3.*p3(3,t3)*A3;
path11 = p3_eq*A3*p2(3,t1)*A2*p3(2,t2)*A3.*p2(3,t3)*A2;
path12 = p3_eq*A3*p2(3,t1)*A2*p3(2,t2)*A3.*p1(3,t3)*A1;
path13 = p3_eq*A3*p2(3,t1)*A2*p2(2,t2)*A2.*p3(2,t3)*A3;
path14 = p3_eq*A3*p2(3,t1)*A2*p2(2,t2)*A2.*p2(2,t3)*A2;
path15 = p3_eq*A3*p2(3,t1)*A2*p2(2,t2)*A2.*p1(2,t3)*A1;
path16 = p3_eq*A3*p2(3,t1)*A2*p1(2,t2)*A1.*p3(1,t3)*A3;
path17 = p3_eq*A3*p2(3,t1)*A2*p1(2,t2)*A1.*p2(1,t3)*A2;
path18 = p3_eq*A3*p2(3,t1)*A2*p1(2,t2)*A1.*p1(1,t3)*A1;

path19 = p3_eq*A3*p1(3,t1)*A1*p3(1,t2)*A3.*p3(3,t3)*A3;
path20 = p3_eq*A3*p1(3,t1)*A1*p3(1,t2)*A3.*p2(3,t3)*A2;
path21 = p3_eq*A3*p1(3,t1)*A1*p3(1,t2)*A3.*p1(3,t3)*A1;
path22 = p3_eq*A3*p1(3,t1)*A1*p2(1,t2)*A2.*p3(2,t3)*A3;
path23 = p3_eq*A3*p1(3,t1)*A1*p2(1,t2)*A2.*p2(2,t3)*A2;
path24 = p3_eq*A3*p1(3,t1)*A1*p2(1,t2)*A2.*p1(2,t3)*A1;
path25 = p3_eq*A3*p1(3,t1)*A1*p1(1,t2)*A1.*p3(1,t3)*A3;
path26 = p3_eq*A3*p1(3,t1)*A1*p1(1,t2)*A1.*p2(1,t3)*A2;
path27 = p3_eq*A3*p1(3,t1)*A1*p1(1,t2)*A1.*p1(1,t3)*A1;

path28 = p2_eq*A2*p3(2,t1)*A3*p3(3,t2)*A3.*p3(3,t3)*A3;
path29 = p2_eq*A2*p3(2,t1)*A3*p3(3,t2)*A3.*p2(3,t3)*A2;
path30 = p2_eq*A2*p3(2,t1)*A3*p3(3,t2)*A3.*p1(3,t3)*A1;
path31 = p2_eq*A2*p3(2,t1)*A3*p2(3,t2)*A2.*p3(2,t3)*A3;
path32 = p2_eq*A2*p3(2,t1)*A3*p2(3,t2)*A2.*p2(2,t3)*A2;
path33 = p2_eq*A2*p3(2,t1)*A3*p2(3,t2)*A2.*p1(2,t3)*A1;
path34 = p2_eq*A2*p3(2,t1)*A3*p1(3,t2)*A1.*p3(1,t3)*A3;
path35 = p2_eq*A2*p3(2,t1)*A3*p1(3,t2)*A1.*p2(1,t3)*A2;
path36 = p2_eq*A2*p3(2,t1)*A3*p1(3,t2)*A1.*p1(1,t3)*A1;

path37 = p2_eq*A2*p2(2,t1)*A2*p3(2,t2)*A3.*p3(3,t3)*A3;
path38 = p2_eq*A2*p2(2,t1)*A2*p3(2,t2)*A3.*p2(3,t3)*A2;
path39 = p2_eq*A2*p2(2,t1)*A2*p3(2,t2)*A3.*p1(3,t3)*A1;
path40 = p2_eq*A2*p2(2,t1)*A2*p2(2,t2)*A2.*p3(2,t3)*A3;
path41 = p2_eq*A2*p2(2,t1)*A2*p2(2,t2)*A2.*p2(2,t3)*A2;
path42 = p2_eq*A2*p2(2,t1)*A2*p2(2,t2)*A2.*p1(2,t3)*A1;
path43 = p2_eq*A2*p2(2,t1)*A2*p1(2,t2)*A1.*p3(1,t3)*A3;
path44 = p2_eq*A2*p2(2,t1)*A2*p1(2,t2)*A1.*p2(1,t3)*A2;
path45 = p2_eq*A2*p2(2,t1)*A2*p1(2,t2)*A1.*p1(1,t3)*A1;

path46 = p2_eq*A2*p1(2,t1)*A1*p3(1,t2)*A3.*p3(3,t3)*A3;
path47 = p2_eq*A2*p1(2,t1)*A1*p3(1,t2)*A3.*p2(3,t3)*A2;
path48 = p2_eq*A2*p1(2,t1)*A1*p3(1,t2)*A3.*p1(3,t3)*A1;
path49 = p2_eq*A2*p1(2,t1)*A1*p2(1,t2)*A2.*p3(2,t3)*A3;
path50 = p2_eq*A2*p1(2,t1)*A1*p2(1,t2)*A2.*p2(2,t3)*A2;
path51 = p2_eq*A2*p1(2,t1)*A1*p2(1,t2)*A2.*p1(2,t3)*A1;
path52 = p2_eq*A2*p1(2,t1)*A1*p1(1,t2)*A1.*p3(1,t3)*A3;
path53 = p2_eq*A2*p1(2,t1)*A1*p1(1,t2)*A1.*p2(1,t3)*A2;
path54 = p2_eq*A2*p1(2,t1)*A1*p1(1,t2)*A1.*p1(1,t3)*A1;

path55 = p1_eq*A1*p3(1,t1)*A3*p3(3,t2)*A3.*p3(3,t3)*A3;
path56 = p1_eq*A1*p3(1,t1)*A3*p3(3,t2)*A3.*p2(3,t3)*A2;
path57 = p1_eq*A1*p3(1,t1)*A3*p3(3,t2)*A3.*p1(3,t3)*A1;
path58 = p1_eq*A1*p3(1,t1)*A3*p2(3,t2)*A2.*p3(2,t3)*A3;
path59 = p1_eq*A1*p3(1,t1)*A3*p2(3,t2)*A2.*p2(2,t3)*A2;
path60 = p1_eq*A1*p3(1,t1)*A3*p2(3,t2)*A2.*p1(2,t3)*A1;
path61 = p1_eq*A1*p3(1,t1)*A3*p1(3,t2)*A1.*p3(1,t3)*A3;
path62 = p1_eq*A1*p3(1,t1)*A3*p1(3,t2)*A1.*p2(1,t3)*A2;
path63 = p1_eq*A1*p3(1,t1)*A3*p1(3,t2)*A1.*p1(1,t3)*A1;

path64 = p1_eq*A1*p2(1,t1)*A2*p3(2,t2)*A3.*p3(3,t3)*A3;
path65 = p1_eq*A1*p2(1,t1)*A2*p3(2,t2)*A3.*p2(3,t3)*A2;
path66 = p1_eq*A1*p2(1,t1)*A2*p3(2,t2)*A3.*p1(3,t3)*A1;
path67 = p1_eq*A1*p2(1,t1)*A2*p2(2,t2)*A2.*p3(2,t3)*A3;
path68 = p1_eq*A1*p2(1,t1)*A2*p2(2,t2)*A2.*p2(2,t3)*A2;
path69 = p1_eq*A1*p2(1,t1)*A2*p2(2,t2)*A2.*p1(2,t3)*A1;
path70 = p1_eq*A1*p2(1,t1)*A2*p1(2,t2)*A1.*p3(1,t3)*A3;
path71 = p1_eq*A1*p2(1,t1)*A2*p1(2,t2)*A1.*p2(1,t3)*A2;
path72 = p1_eq*A1*p2(1,t1)*A2*p1(2,t2)*A1.*p1(1,t3)*A1;

path73 = p1_eq*A1*p1(1,t1)*A1*p3(1,t2)*A3.*p3(3,t3)*A3;
path74 = p1_eq*A1*p1(1,t1)*A1*p3(1,t2)*A3.*p2(3,t3)*A2;
path75 = p1_eq*A1*p1(1,t1)*A1*p3(1,t2)*A3.*p1(3,t3)*A1;
path76 = p1_eq*A1*p1(1,t1)*A1*p2(1,t2)*A2.*p3(2,t3)*A3;
path77 = p1_eq*A1*p1(1,t1)*A1*p2(1,t2)*A2.*p2(2,t3)*A2;
path78 = p1_eq*A1*p1(1,t1)*A1*p2(1,t2)*A2.*p1(2,t3)*A1;
path79 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1.*p3(1,t3)*A3;
path80 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1.*p2(1,t3)*A2;
path81 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1.*p1(1,t3)*A1;
% disp(['Size path81 = ' num2str(size(path81))]);

kappa2 = path1 + path2 + path3 + path4 + path5 + path6 + path7 + path8 + path9 + path10 + path11 + path12 + path13 + path14 + path15 + path16 + path17 + path18 + path19 + path20;
kappa2 = kappa2 + path21 + path22 + path23 + path24 + path25 + path26 + path27 + path28 + path29 + path30 + path31 + path32 + path33 + path34 + path35 + path36 + path37 + path38 + path39 + path40;
kappa2 = kappa2 + path41 + path42 + path43 + path44 + path45 + path46 + path47 + path48 + path49 + path50 + path51 + path52 + path53 + path54 + path55 + path56 + path57 + path58 + path59 + path60;
C4 = kappa2 + path61 + path62 + path63 + path64 + path65 + path66 + path67 + path68 + path69 + path70 + path71 + path72 + path73 + path74 + path75 + path76 + path77 + path78 + path79 + path80 + path81;
% disp(['Size C4 = ' num2str(size(C4))]);

C2_left = (p3_eq*A3*(A3*p3(3,t1) + A2*p2(3,t1) + A1*p1(3,t1)) + p2_eq*A2*(A3*p3(2,t1) + A2*p2(2,t1) + A1*p1(2,t1)) + p1_eq*A1*(A3*p3(1,t1) + A2*p2(1,t1) + A1*p1(1,t1)));
C2 = C2_left;
C2_right = (p3_eq*A3*(A3*p3(3,t3) + A2*p2(3,t3) + A1*p1(3,t3)) + p2_eq*A2*(A3*p3(2,t3) + A2*p2(2,t3) + A1*p1(2,t3)) + p1_eq*A1*(A3*p3(1,t3) + A2*p2(1,t3) + A1*p1(1,t3)));
C2Product = C2_left.*C2_right;
% disp(['Size C2Product = ' num2str(size(C2Product))]);
C4_diff = C4 - C2Product;

% disp(['Size C4_diff = ' num2str(size(C4_diff))]);

 %-------------------------------------------------------------------------
 % Plot the data
 %-------------------------------------------------------------------------
 if plotMode == 1
        clf;
        surf(t1, t3, C4);
        title('Four-point TCF: C^{(4)}','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('Time (\tau_3)','FontSize',14);
        zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
        
        view(28,36);
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
        
        drawnow();
        hold on;
 end
    
 %-------------------------------------------------------------------------
 %
 %-------------------------------------------------------------------------
    function p1 = p1(init,time)
        if init == 3
            m = m_3;
            n = n_3;
        elseif init == 2
            m = m_2;
            n = n_2;
        else
            m = m_1;
            n = n_1;
        end
        
        p1 = p1_eq + m*a*exp(-lam1*time) + n*x*exp(-lam2*time);
    end

    function p2 = p2(init,time)
        if init == 3
            m = m_3;
            n = n_3;
        elseif init == 2
            m = m_2;
            n = n_2;
        else
            m = m_1;
            n = n_1;
        end
        
        p2 = p2_eq + m*b*exp(-lam1*time) + n*y*exp(-lam2*time);
    end

    function p3 = p3(init,time)
        if init == 3
            m = m_3;
            n = n_3;
        elseif init == 2
            m = m_2;
            n = n_2;
        else
            m = m_1;
            n = n_1;
        end
        
        p3 = p3_eq + m*c*exp(-lam1*time) + n*z*exp(-lam2*time);
    end
end


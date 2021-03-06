function C2_sim = C2maker_3state123_cyclical_analytical(t12,t13,t21,t23,t31,A1,A2,A3,time)
% Now is a mean subtracted 2-point TCF. Works the same as TCF_cyclic3state
programName = 'C2maker_3state123_cyclical_analytical';
%--------------------------------------------------------------------------
% SET PARAMATERS
%--------------------------------------------------------------------------
switch nargin
    case 0
        disp(['Using default values in ' programName]);
        
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
        
        Nparams = length(boundsArray);
        population = rand(1,Nparams);
        for param_idx = 1:Nparams
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
        
        Npts = 150;
        time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
    case 8
        
        Npts = 150;
        time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
        
end

% note: this messes it up:
%t12 = 0.000825, t13 = 0.006767, t21 = 0.000615, t23 = 0.006725, t31 = 0.003205, t32 = 0.004270  A1 = 0.737326, A2 = 0.518160, A3 = 0.315085,
%--------------------------------------------------------------------------
% Set the rates
%--------------------------------------------------------------------------
k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;
% Detailed balance condition:
k32 = k12*k23*k31/(k13*k21);

% % OUTOUT FROM ODEsolver_3state123_cyclical
% % Eigenvalue 1 = 0
% % Eigenvalue 2 = - k12/2 - k13/2 - k21/2 - k23/2 - k31/2 - k32/2 - (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2
% % Eigenvalue 3 = (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2 - k13/2 - k21/2 - k23/2 - k31/2 - k32/2 - k12/2

% Define eigenvalues lam1 and lam2 (3state123_cyclical_analytical)
c1 = 0.5*(k12+k21+k23+k32+k13+k31);
% c2 = 0.5*sqrt(k12^2 + 2*k12*k13 + 2*k12*k21 - 2*k12*k23 - 2*k12*k31 - 2*k12*k32 + k13^2 - 2*k13*k21 - 2*k13*k23 + 2*k13*k31 - k13*k32 + k21^2 + 2*k21*k23 - 2*k21*k31 - 2*k21*k32 + k23^2 - 2*k23*k31 + 2*k23*k32 + k31^2 + 2*k31*k32 + k32^2);
c2 = 0.5*sqrt(k12^2 + 2*k12*k13 + 2*k12*k21 - 2*k12*k23 - 2*k12*k31 - 2*k12*k32 + k13^2 - 2*k13*k21 - 2*k13*k23 + 2*k13*k31 - 2*k13*k32 + k21^2 + 2*k21*k23 - 2*k21*k31 - 2*k21*k32 + k23^2 - 2*k23*k31 + 2*k23*k32 + k31^2 + 2*k31*k32 + k32^2);

lam1 = c1 + c2;
lam2 = c1 - c2;
% disp(['1/Lam 1 = ' num2str(1/lam1) ]);
% disp(['1/Lam 2 = ' num2str(1/lam2) ]);

% % OUTOUT FROM ODEsolver_3state123_cyclical: Same values
% % Eigenvector2_1 = (k23 + k31 + k32)/(k13 - k23) - (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 + (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23)
% % Eigenvector2_2 = (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 + (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23) - (k13 + k31 + k32)/(k13 - k23)
% % Eigenvector2_3 = 1
% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
a = (k23 + k31 + k32 - c1 - c2)/(k13 - k23);
b = (c1 + c2 - k13 - k31 - k32)/(k13 - k23);
c = 1;

% % OUTOUT FROM ODEsolver_3state123_cyclical: Same values
% % Eigenvector3_1 = (k23 + k31 + k32)/(k13 - k23) - (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 - (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23)
% % Eigenvector3_2 = (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 - (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23) - (k13 + k31 + k32)/(k13 - k23)
% % Eigenvector3_3 = 1
% Define eigenvector components corresponding to lam2: v2 = [x,y,z].
x = (k23 + k31 + k32 - c1 + c2)/(k13 - k23);
y = (c1 - c2 - k13 - k31 - k32)/(k13 - k23);
z = 1;

% % Using PeqSolver_3state123_cyclical (solve): we get same ansers
% P1eqSol = (k21*k31 + k21*k32 + k23*k31)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)
%  P2eqSol = (k12*k31 + k12*k32 + k13*k32)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)
% P3eqSol = (k13*k21 + k12*k23 + k13*k23)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)

% Define equilibrium values of populations
p1_eq = 1/(1 + (k12 + k13)/k21 + (1 - k31/k21)*((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21));
p3_eq = ((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21)*p1_eq;
p2_eq = ((k12 + k13)*p1_eq - k31*p3_eq)/k21;

% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_3 = (-p2_eq + b/a*p1_eq)/(y - b/a*x);
n_2 = (1 - p2_eq + b/a*p1_eq)/(y - b/a*x);
n_1 = (-p2_eq - b/a*(1 - p1_eq))/(y - b/a*x);
m_3 = (-p1_eq - n_3*x)/a;
m_2 = (-p1_eq - n_2*x)/a;
m_1 = (1 - p1_eq - n_1*x)/a;


%%Using ODE solver: GIVES THE SAME ANSER!
% t = time;
% tic
% p33 = (k13*k21 + k12*k23 + k13*k23)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31) + (0.5*exp(-0.5*t*(k12 + k13 + k21 + k23 + k31 + k32 + (2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)))*(k12*k31^2 - 1.0*k12^2*k31 + k12*k32^2 - 1.0*k12^2*k32 + k13*k32^2 - 1.0*k13^2*k32 + k21*k31^2 - 1.0*k21^2*k31 + k21*k32^2 - 1.0*k21^2*k32 + k23*k31^2 - 1.0*k23^2*k31 + k12*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k12*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k13*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k23*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) - 1.0*k12*k13*k31 - 2.0*k12*k13*k32 - 2.0*k12*k21*k31 - 2.0*k12*k21*k32 + k13*k21*k31 + k12*k23*k32 + k13*k23*k31 + k13*k23*k32 + 2.0*k12*k31*k32 - 2.0*k21*k23*k31 + k13*k31*k32 - 1.0*k21*k23*k32 + 2.0*k21*k31*k32 + k23*k31*k32))/((k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)) + (0.5*exp(-0.5*t*(k12 + k13 + k21 + k23 + k31 + k32 - 1.0*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)))*(k12^2*k31 - 1.0*k12*k31^2 - 1.0*k12*k32^2 + k12^2*k32 - 1.0*k13*k32^2 + k13^2*k32 - 1.0*k21*k31^2 + k21^2*k31 - 1.0*k21*k32^2 + k21^2*k32 - 1.0*k23*k31^2 + k23^2*k31 + k12*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k12*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k13*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k23*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k12*k13*k31 + 2.0*k12*k13*k32 + 2.0*k12*k21*k31 + 2.0*k12*k21*k32 - 1.0*k13*k21*k31 - 1.0*k12*k23*k32 - 1.0*k13*k23*k31 - 1.0*k13*k23*k32 - 2.0*k12*k31*k32 + 2.0*k21*k23*k31 - 1.0*k13*k31*k32 + k21*k23*k32 - 2.0*k21*k31*k32 - 1.0*k23*k31*k32))/((k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2));
% toc

% Define populations as a function of time for different initial values.
% p2_1 population 2 as a function of time, with the population initially in
% state 1.
% tic
p3_3 = p3_eq + m_3*c*exp(-lam1*time) + n_3*z*exp(-lam2*time);
% toc
p3_2 = p3_eq + m_2*c*exp(-lam1*time) + n_2*z*exp(-lam2*time);
p3_1 = p3_eq + m_1*c*exp(-lam1*time) + n_1*z*exp(-lam2*time);
p2_3 = p2_eq + m_3*b*exp(-lam1*time) + n_3*y*exp(-lam2*time);
p2_2 = p2_eq + m_2*b*exp(-lam1*time) + n_2*y*exp(-lam2*time);
p2_1 = p2_eq + m_1*b*exp(-lam1*time) + n_1*y*exp(-lam2*time);
p1_3 = p1_eq + m_3*a*exp(-lam1*time) + n_3*x*exp(-lam2*time);
p1_2 = p1_eq + m_2*a*exp(-lam1*time) + n_2*x*exp(-lam2*time);
p1_1 = p1_eq + m_1*a*exp(-lam1*time) + n_1*x*exp(-lam2*time);

% Subtract mean values
Amean = p1_eq*A1 + p2_eq*A2 + p3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

% Calculate C2_sim
C2_sim = p3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
    p2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
    p1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
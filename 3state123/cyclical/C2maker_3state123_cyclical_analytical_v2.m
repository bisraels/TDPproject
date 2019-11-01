% function C2_sim = C2maker_3state123_cyclical_analytical_v2(t12,t13,t21,t23,t31,A1,A2,A3,time)
% Now is a mean subtracted 2-point TCF. Works the same as TCF_cyclic3state
programName = 'C2maker_3state123_cyclical_analytical_v2';
plotMode = 1;
%--------------------------------------------------------------------------
% SET PARAMATERS
%--------------------------------------------------------------------------
% switch nargin
%     case 0
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
        for parac1_idx = 1:Nparams
            %To pick a random number in the interval of LB to UB:
            % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
            population(parac1_idx) = boundsArray(parac1_idx) + population(parac1_idx)*(boundsArray(parac1_idx,2) - boundsArray(parac1_idx,1));
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
%     case 8
%         
%         Npts = 150;
%         time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
%         
% end

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

% % Using PeqSolver_3state123_cyclical (solve): we get same ansers
% p1_eq = (k21*k31 + k21*k32 + k23*k31)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)
%  p2_eq = (k12*k31 + k12*k32 + k13*k32)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)
% p3_eq = (k13*k21 + k12*k23 + k13*k23)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)

% Define equilibrium values of populations
P1_eq = 1/(1 + (k12 + k13)/k21 + (1 - k31/k21)*((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21));
P3_eq = ((k21 + k23)*(k12 + k13) - k12*k21)/((k21 + k23)*k31 + k32*k21)*P1_eq;
P2_eq = ((k12 + k13)*P1_eq - k31*P3_eq)/k21;

% % OUTOUT FROM ODEsolver_3state123_cyclical:
% % Eigenvector1_1 = (k21*k31 + k21*k32 + k23*k31)/(k13*k21 + k12*k23 + k13*k23)
% % Eigenvector1_2 = (k12*k31 + k12*k32 + k13*k32)/(k13*k21 + k12*k23 + k13*k23)
% % Eigenvector1_3 = 1

% % OUTOUT FROM ODEsolver_3state123_cyclical: Same values
% % Eigenvector2_1 = (k23 + k31 + k32)/(k13 - k23) - (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 + (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23)
% % Eigenvector2_2 = (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 + (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23) - (k13 + k31 + k32)/(k13 - k23)
% % Eigenvector2_3 = 1
% Define eigenvector components corresponding to lam2: v2 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
v2_1 = (k23 + k31 + k32 - c1 - c2)/(k13 - k23);
v2_2 = (c1 + c2 - k13 - k31 - k32)/(k13 - k23);
v2_3 = 1;

% % OUTOUT FROM ODEsolver_3state123_cyclical: Same values
% % Eigenvector3_1 = (k23 + k31 + k32)/(k13 - k23) - (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 - (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23)
% % Eigenvector3_2 = (k12/2 + k13/2 + k21/2 + k23/2 + k31/2 + k32/2 - (2*k12*k13 + 2*k12*k21 - 2*k13*k21 - 2*k12*k23 - 2*k13*k23 - 2*k12*k31 - 2*k12*k32 + 2*k13*k31 + 2*k21*k23 - 2*k13*k32 - 2*k21*k31 - 2*k21*k32 - 2*k23*k31 + 2*k23*k32 + 2*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)/2)/(k13 - k23) - (k13 + k31 + k32)/(k13 - k23)
% % Eigenvector3_3 = 1
% Define eigenvector components corresponding to lam3: v2 = [x,y,z].
v3_1 = (k23 + k31 + k32 - c1 + c2)/(k13 - k23);
v3_2 = (c1 - c2 - k13 - k31 - k32)/(k13 - k23);
v3_3 = 1;


% Define constants needed for differential equation solutions, for various
% initial values of the population. c1_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
c3_3 = (-P2_eq + v2_2/v2_1*P1_eq)/(v3_2 - v2_2/v2_1*v3_1);
c3_2 = (1 - P2_eq + v2_2/v2_1*P1_eq)/(v3_2 - v2_2/v2_1*v3_1);
c3_1 = (-P2_eq - v2_2/v2_1*(1 - P1_eq))/(v3_2 - v2_2/v2_1*v3_1);
c2_3 = (-P1_eq - c3_3*v3_1)/v2_1;
c2_2 = (-P1_eq - c3_2*v3_1)/v2_1;
c2_1 = (1 - P1_eq - c3_1*v3_1)/v2_1;


%%Using ODE solver: GIVES THE SAME ANSER!
% t = time;
% tic
% p33 = (k13*k21 + k12*k23 + k13*k23)/(k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31) + (0.5*exp(-0.5*t*(k12 + k13 + k21 + k23 + k31 + k32 + (2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)))*(k12*k31^2 - 1.0*k12^2*k31 + k12*k32^2 - 1.0*k12^2*k32 + k13*k32^2 - 1.0*k13^2*k32 + k21*k31^2 - 1.0*k21^2*k31 + k21*k32^2 - 1.0*k21^2*k32 + k23*k31^2 - 1.0*k23^2*k31 + k12*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k12*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k13*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k23*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) - 1.0*k12*k13*k31 - 2.0*k12*k13*k32 - 2.0*k12*k21*k31 - 2.0*k12*k21*k32 + k13*k21*k31 + k12*k23*k32 + k13*k23*k31 + k13*k23*k32 + 2.0*k12*k31*k32 - 2.0*k21*k23*k31 + k13*k31*k32 - 1.0*k21*k23*k32 + 2.0*k21*k31*k32 + k23*k31*k32))/((k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)) + (0.5*exp(-0.5*t*(k12 + k13 + k21 + k23 + k31 + k32 - 1.0*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2)))*(k12^2*k31 - 1.0*k12*k31^2 - 1.0*k12*k32^2 + k12^2*k32 - 1.0*k13*k32^2 + k13^2*k32 - 1.0*k21*k31^2 + k21^2*k31 - 1.0*k21*k32^2 + k21^2*k32 - 1.0*k23*k31^2 + k23^2*k31 + k12*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k12*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k13*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k21*k32*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k23*k31*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2) + k12*k13*k31 + 2.0*k12*k13*k32 + 2.0*k12*k21*k31 + 2.0*k12*k21*k32 - 1.0*k13*k21*k31 - 1.0*k12*k23*k32 - 1.0*k13*k23*k31 - 1.0*k13*k23*k32 - 2.0*k12*k31*k32 + 2.0*k21*k23*k31 - 1.0*k13*k31*k32 + k21*k23*k32 - 2.0*k21*k31*k32 - 1.0*k23*k31*k32))/((k13*k21 + k12*k23 + k13*k23 + k12*k31 + k12*k32 + k13*k32 + k21*k31 + k21*k32 + k23*k31)*(2.0*k12*k13 + 2.0*k12*k21 - 2.0*k13*k21 - 2.0*k12*k23 - 2.0*k13*k23 - 2.0*k12*k31 - 2.0*k12*k32 + 2.0*k13*k31 + 2.0*k21*k23 - 2.0*k13*k32 - 2.0*k21*k31 - 2.0*k21*k32 - 2.0*k23*k31 + 2.0*k23*k32 + 2.0*k31*k32 + k12^2 + k13^2 + k21^2 + k23^2 + k31^2 + k32^2)^(1/2));
% toc

% Define populations as a function of time for different initial values.
% p2_1 population 2 as a function of time, with the population initially in
% state 1. 
% tic
% toc
p1_1 = P1_eq + c2_1*v2_1*exp(-lam1*time) + c3_1*v3_1*exp(-lam2*time);
p1_2 = P1_eq + c2_2*v2_1*exp(-lam1*time) + c3_2*v3_1*exp(-lam2*time);
p1_3 = P1_eq + c2_3*v2_1*exp(-lam1*time) + c3_3*v3_1*exp(-lam2*time);
p2_1 = P2_eq + c2_1*v2_2*exp(-lam1*time) + c3_1*v3_2*exp(-lam2*time);
p2_2 = P2_eq + c2_2*v2_2*exp(-lam1*time) + c3_2*v3_2*exp(-lam2*time);
p2_3 = P2_eq + c2_3*v2_2*exp(-lam1*time) + c3_3*v3_2*exp(-lam2*time);
p3_1 = P3_eq + c2_1*v2_3*exp(-lam1*time) + c3_1*v3_3*exp(-lam2*time);
p3_2 = P3_eq + c2_2*v2_3*exp(-lam1*time) + c3_2*v3_3*exp(-lam2*time);
p3_3 = P3_eq + c2_3*v2_3*exp(-lam1*time) + c3_3*v3_3*exp(-lam2*time);

% Subtract mean values
Amean = P1_eq*A1 + P2_eq*A2 + P3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

% Calculate C2_sim
tic
C2_sim = P3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
    P2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
    P1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
toc
% Matrix of conditional probabilities Pji (i-->j) with i is initial condition
% cP = [ p1_1, p1_2, p1_3;
%     p2_1, p2_2, p2_3;
%     p3_1, p3_2, p3_3];

%mean has already been subtracted
A = [A1,A2,A3];
Peq = [P1_eq,P2_eq,P3_eq];

cP = zeros(numel(A),numel(A),length(time));
cP(1,1,:) = p1_1;
cP(2,1,:) = p2_1;
cP(3,1,:) = p3_1;
cP(1,2,:) = p1_2;
cP(2,2,:) = p2_2;
cP(3,2,:) = p3_2;
cP(1,3,:) = p1_3;
cP(2,3,:) = p2_3;
cP(3,3,:) = p3_3;

tic
C2_sim2 = zeros(size(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        C2_sim1_temp = A(j) * squeeze(cP(j,i,:)) * A(i) * Peq(i);
        C2_sim2 = C2_sim2 + C2_sim1_temp;
    end
end
disp('Time to calculate C2 using loops...');
toc
%--------------------------------------------------------------------------
%  Plot two point TCF
%--------------------------------------------------------------------------

%close all
if plotMode == 1
    figure(2)
    clf;
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    %subplot(1,2,1)
    %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
    TCF2pt = plot(time,C2_sim,'LineWidth',2);
    hold on;
     TCF2pt2 = plot(time,C2_sim2,'r--','LineWidth',2);
    
    title('Two point TCF','FontSize',18)
    xlabel('Time','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);
    %xlim([10^-5 10^0])
    % ylim([C2(500) C2(0)])
    
    ax = gca;
    ax.XScale = 'log';
%     set(gca,'yscale','log')
    
end


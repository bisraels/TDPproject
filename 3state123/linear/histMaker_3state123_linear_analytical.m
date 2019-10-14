function [Peq] = histMaker_3state123_linear_analytical(t12,t21,t23,t32,A1,A2,A3)
%MODEL: 1 <--> 2 <--> 3
programName = 'histMaker_3state123_linear_analytical.m';
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
end

%--------------------------------------------------------------------------
% Set the rates (4rates)
%--------------------------------------------------------------------------
k12 = 1/t12;
k21 = 1/t21;
k23 = 1/t23;
k32 = 1/t32;

% Define eigenvalues lam1 and lam2\  
% c1 = 0.5*(k12+k21+k23+k32);
% c2 = 0.5*sqrt(k12^2 + 2*k12*k21 - 2*k12*k23 - 2*k12*k32 + k21^2 + 2*k21*k23 - 2*k21*k32 + k23^2 + 2*k23*k32 + k32^2);
% lam1 = c1 + c2;
% lam2 = c1 - c2;

% Define equilibrium values of populations
p2_eq = 1/(k21/k12 + 1 + k23/k32);
p1_eq = k21/k12*p2_eq;
p3_eq = k23/k32*p2_eq;

% Define equilibrium values of populations
Peq = [p1_eq; p2_eq; p3_eq];

end


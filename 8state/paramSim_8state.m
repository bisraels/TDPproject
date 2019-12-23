function [K,A,time,rates] = paramSim_8state()
%  Need to set bounds for each time constant
% Set DNA-only states from data fitting of 3state model - 3p15mer
t12_bounds = [65e-6,85e-6];        %Paramater #1 is high--> med
t13_bounds = [16e-3,18e-3];       %Paramater #2 is high --> low
t21_bounds = [25e-6,29e-6];           %Paramater #3 is med --> high
t23_bounds = [1e-3,200e-3];%***          %Paramater #4 is med --> low
t31_bounds = [2.5e-4,3.5e-4];         %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
% *t32 wll be determined by the others in loop

% Protein binding time constants:
t24_bounds = [100e-6, 1e-1];
t42_bounds = [10e-6, 1e-1];
t25_bounds = [10e-6, 1e-1];
t52_bounds = [10e-6, 1e-1];
t56_bounds = [10e-6, 1e-1];
% * t65_bounds will be determined by others in loop
t63_bounds = [10e-6, 1e-1];
t36_bounds = [10e-6, 1e-1];
t37_bounds = [10e-6, 1e-1];
t73_bounds = [10e-6, 1e-1];
t68_bounds = [10e-6, 1e-1];
t86_bounds = [10e-6, 1e-1];
% t89_bounds = [1e-6, 1e-1];
% t98_bounds = [1e-6, 1e-1];
% % * t96_bounds will be determined by others in loop
% t69_bounds = [1e-6, 1e-1];

% Define bounds of FRET parameters
A1_bounds = [0.765,0.775];   % Compact       (higest FRET)            %Paramater #6 % HIGH fret State
A2_bounds = [0.54,0.56]; % Not extended  (Intermediate FRET)      %Paramater #7 % Med FRET state
A3_bounds = [0.1,0.25];    % Extended      (lowest FRET)            %Paramater #8 %Low FRET state
A4_bounds = [0.3,0.65]; % Not extended
A5_bounds = [0.3,0.65]; % Not extended
A6_bounds = [0,0.25];    % Extended
A7_bounds = [0,0.25];    % Extended
A8_bounds = [0,0.25];    % Extended
% A9_bounds = [0,0.3]


boundsArray = [t12_bounds; t13_bounds; t21_bounds; t23_bounds; t31_bounds;...
    t24_bounds; t42_bounds; t25_bounds; t52_bounds; t56_bounds; t63_bounds;...
    t36_bounds; t37_bounds; t73_bounds; t68_bounds; t86_bounds;...
    A1_bounds; A2_bounds; A3_bounds; A4_bounds; A5_bounds;A6_bounds; A7_bounds;...
    A8_bounds;];

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
% * t32_bounds will be determined by others in loop
t24 = population(6);
t42 = population(7);
t25 = population(8);
t52 = population(9);
t56 = population(10);
t65 = population(11);
% * t63_bounds will be determined by others in loop
t36 = population(12);
t37 = population(13);
t73 = population(14);
t68 = population(15);
t86 = population(16);
%             t89 = population(pop_idx,17);
%             t98 = population(pop_idx,18);
% * t96_bounds will be determined by others in loop
%             t69 = population(pop_idx,20);

%             % Loop conditions
%             t32 = (t12 * t23 * t31)/(t13 * t21);
%             t63 = (t23 * t36 * t65* t52)/(t32 * t25 * t56);
%
% Need these indexes for 9 state
A1 = population(17);        %21);
A2 = population(18);        %22);
A3 = population(19);        %23);
A4 = population(20);        %24);
A5 = population(21);        %25);
A6 = population(22);        %26);
A7 = population(23);        %27);
A8 = population(24);        %28);
%             A9 = population(pop_idx,29);

% Define rates, kij, from the tij's
k12 = 1/t12;
k13 = 1/t13;
k21 = 1/t21;
k23 = 1/t23;
k31 = 1/t31;
k24 = 1/t24;
k42 = 1/t42;
k25 = 1/t25;
k52 = 1/t52;
k56 = 1/t56;
k65 = 1/t65;
k36 = 1/t36;
k37 = 1/t37;
k73 = 1/t73;
k68 = 1/t68;
k86 = 1/t86;
%           k89 = 1/t89;
%           k98 = 1/t98;
%           t96 = 1/t96;
%           t69 = 1/t69;

% Loop conditions
k32 = (k12 * k23 * k31)/(k13 * k21);
k63 = (k23 * k36 * k65* k52)/(k32 * k25 * k56);
t32 = 1/k32;
t63 = 1/k63;

%             A = [A1;A2;A3;A4;A5;A6;A7;A8];
A = [A1,A2,A3,A4,A5,A6,A7,A8];
rates = [k12,k13,k21,k23,k32,k31,k24,k42,k25,k52,k56,k65,k63,k36,k37,k73,k68,k86];

% Define K matrix
% 8 state model: VERIFIED CORRECT
K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0;...
    k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0;...
    k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0;...
    0, k24, 0, -k42, 0, 0, 0, 0;...
    0, k25, 0, 0, -(k52+k56), k65, 0, 0;...
    0, 0, k36, 0, k56, -(k65+k63+k68), 0, k86;...
    0, 0, k37, 0, 0, 0, -k73, 0;...
    0, 0, 0, 0, 0, k68, 0, -k86;];

Npts = 150;
time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;


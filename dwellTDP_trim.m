close all
clear all
clc

%%
% Claire Albrecht
% July 2019

% Takes dwell functions of HMM fits (made in pathdwellredo) and turn it into TDPs.
% Using Brett's file (dwell2tdp_fromBrett) as a template to understand how the program works

% Dwell file format:
% FRET_initial      FRET_final     dwelltime (in frames)
% to get time from dwell time, multiply by resolution.

% NOTE: This is the second version of this file, I got lost in the weeds
% the first time and wanted to start fresh.


% Set FRET limits
FRET_LB = 0;
FRET_UB = 1;

% Find dwell files
scanID = '*';
fileKeyWord = [scanID 'dwell.dat'];
fileNames = dir(fileKeyWord);
numFiles = length(fileNames);

% Define parameters
tau = 30000;                    % Vaule of tau in ms - figure out how to decide this later.
sigmaFRET = 0.015;              % Width of gaussian
sigma_str = num2str(sigmaFRET); % Converts numbers to words
res = 0.03;               % Later go back and define resolution from data file name - more efficient.

% Define parameters for the histogram
gridRes = 1001;  % Number of bins to organize the grid into
gridMin = 0;    % Minimum value of grid - needs to be 0
gridMax = 1;    % Maximum value of grid - needs to be 1

% Using these parameters define dimensions of grid [min:stepvalue:max] - check this
gridDimension = gridRes.*[gridMin:(gridMax-gridMin)/(gridRes-1):gridMax];


%%


% Initialize matrix to hold 2D transition matrix values
tdpGrid = zeros(gridRes);   % matrix of zeros with gridRes = n for nxn matrix
% Define FRET arrays: fx (fretInitial), fy(fretFinal)
fx = 0:0.001:1; % initial:step value:final
fy = fx;
% Define X & Y as arrays (matricies) that will allow us to operate fx
% and fy like a functions
[X,Y] = meshgrid(fx,fy); % meshgrid gives a 2D grid of coordinates from coordinates in fx and fy.

% Number of transitions recorded. Reset to zero before opening files.
N_transitions_total = 0;
% Number of transitions that happen faster than tau
N_legal_transitions = 0;
N_legal_files = 0;

% Look at molecules one at a time until the end of the list
for molID = 1:numFiles
    color_temp = [rand rand rand];
    dwellFileName = fileNames(molID).name;
    
    % disp(molID)
    
    %disp('for loop to choose molecule file')

    % scanID is the unique identifier of the trace. Take it from  the
    % file names (return to this later).
    scanIDFormat = 'YYYYMMDD-001';
    scanID = dwellFileName(1:length(scanIDFormat));
    
    
    if fileNames(molID).bytes > 0
        ind_tdpGrid = zeros(gridRes);            % Create grid of zeros with dimensions of gridRes.
        N_legal_files = N_legal_files + 1;       % redefine N_legal_files
        dwelltimes = dlmread(dwellFileName);     % Call dwelltimes from the dwell files

        FRET_initial = dwelltimes(:,1);          % Take data from the proper columns in the dwell files
        FRET_final = dwelltimes(:,2);
        dwell_time = dwelltimes(:,3);
        
        %disp(length(dwell_time));                % Show me the length of the run for each molecule
        %disp('for');
        %disp(molID);                             % Show me for which molecule
        
        %disp('counting allowed files')

        hmmFit = [];                             % Empty array
        [Ntransitions,~] = size(dwelltimes);     % Rows containing number of transitions
        timept = 1;
        
        for n = 1:Ntransitions
            N_transitions_total = N_transitions_total + 1;
            hmm_Fit(timept:timept + dwell_time(n)) = FRET_initial(n);
            timept = timept + dwell_time(n);
            
            % Check if transition was shorter than the cut off
            if dwell_time(n) * res*10^-2 <= tau
                N_legal_transitions = N_legal_transitions + 1; % If transition is less than cutoff, count it.
                fret_initial_temp = FRET_initial(n) * gridRes;
                fret_final_temp = FRET_final(n) * gridRes;
                tempIndex_fret_initial_temp = find(gridDimension <= fret_initial_temp,1,'last');
                tempIndex_fret_final_temp = find(gridDimension <= fret_final_temp,1,'last'); % find locates positions where condition is satisfied? dont totally understand this)
                
                % Add a transition to total TDP Grid
                tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) = tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) + 1;
                ind_tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) = ind_tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) + 1;
                % Brett has other stuff I am going to leave for now
                
                %disp('check if transition is shorter than cutoff')
                
            end
            
            % Want to eliminate transition if it occurs at the end of one
            % file, where it is stitched to another.
            if n == Ntransitions
                N_transitions_total = N_transitions_total - 1;
                
                tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) = tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) - 1;
                ind_tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) = ind_tdpGrid(tempIndex_fret_initial_temp,tempIndex_fret_final_temp) - 1;
                
            end

        end
    end
end

%%
%  Make TDP for all molecules

tdpGrid_gauss = zeros(size(tdpGrid)); % An nxn array of zeros where n = size of tdpGrid

index_array = find(tdpGrid > 0);    % the find function gives the positions in the array where the condition is satisfied
Npts = length(index_array);         % This will tell us how many of the elements of the array satisfy the condition
% s is the size of the array

if sigmaFRET > 0
    for i=1:Npts
        [x0_index,y0_index] = ind2sub(size(tdpGrid),index_array(i));
        x0 = fx(x0_index);
        y0 = fy(y0_index);
        A = tdpGrid(index_array(i)); % Amplitude of gaussian
        % Generate a gaussian centered at (x0, y0) with amplitude A.
        gaussianPeak = A * exp(-( (X - x0).^2 + (Y - y0).^2 ) / (2 * sigmaFRET^2) );
        % Add the gaussian at point i to the TDP
        tdpGrid_gauss = tdpGrid_gauss + gaussianPeak;
        
        
    end
    
else
    % Gaussian weighted TDP is equal to the TDP if sigma = 0;
    tdpGrid_gauss = tdoGrid;
end

%  Plot a single  TDP

figure
tdpGrid_gauss_normalized = tdpGrid_gauss./sum(tdpGrid_gauss(:));
hold on
fig1 = surf(X,Y,tdpGrid_gauss_normalized,'FaceAlpha',0.5);
%fig1 = surf(X,Y,tdpGrid_gauss_normalized)
%xlim and ylim
 shading interp % this emliminates teh black lines from the gridding

%%

% Save figure
foutName = [ 'TDP_' scanID 'res' num2str(res) 's_' 'tau' num2str(tau) 'ms'];
extension = 'png';
folderOutName = 'dwellTDP_2_output';
if exist(folderOutName,'dir') ~= 7
    mkdir(folderOutName);
    disp(['Creating a folder named' folderOutName 'to store output']);
end


%saveas(fig1,[folderOutName filesep() foutName],'fig');

saveas(fig1, 'TDP_3p15mer')









%% Find local maxima and values of peaks

L = length(tdpGrid);
Z = tdpGrid_gauss_normalized;   % Z-axis of the surface plot

Num_loc_max = 0;

%mtx_trns = zeros(1,3);

for i = 2:(L-1)
    for j = 2:(L-1)
        
        cell = Z(i,j);          % Zis the matrix that we are plotting in the surface plot.
        
        % Evaluate neighboring matrix cells
        above = Z(i, j + 1);
        below = Z(i, j - 1);
        left = Z(i - 1, j);
        right = Z(i + 1,j);
        
        if above < cell && below < cell && left < cell && right < cell
            Num_loc_max = Num_loc_max + 1;
           % disp('Local max found at (x ,y, z) =');
           % disp([(i)/201,(j)/201, Z(i,j)]);     
            
            if exist('mtx_trns')
                mtx_trns = vertcat(mtx_trns, [i/gridRes,j/gridRes, Z(i,j)]);
            else
                mtx_trns = [(i)/gridRes,(j)/gridRes, Z(i,j)];
            end
        end
    end
end

Num_loc_max;
mtx_trns;

%% Try plotting FRET values from HMM - THIS DID NOT WORK.
% But it is helpful to have hmmFRET defined in this file.

fileEndName = 'PATH.dat';
fileNames = dir(['*' fileEndName]);
numFiles = length(fileNames);

for molNum = 1:numFiles  
    fileName = fileNames(molNum).name;   
    pathFile = importdata(fileName,' ',0);   
    hmmFRET = pathFile(:,5);
end


FRET_vals = unique(hmmFRET);
plot3(mtx_trns(:,1),mtx_trns(:,2),mtx_trns(:,3),'r.','MarkerSize',14)
% hold on
% plot(FRET_vals(1),FRET_vals(2),1.8*10^-3,'blue.','MarkerSize',15)
% hold on
% plot(FRET_vals(2),FRET_vals(3),1.8*10^-3,'blue.','MarkerSize',15)
% hold on
% plot(FRET_vals(3),FRET_vals(4),1.8*10^-3,'blue.','MarkerSize',15)
% hold on
% plot(FRET_vals(4),FRET_vals(5),1.8*10^-3,'blue.','MarkerSize',15)
hold off



%% Identify the rates that exist based on the TDP

syms k12 k13 k14 k15 k21 k23 k24 k25 k31 k32 k34 k35 k41 k42 k43 k45 k51 k52 k53 k54 


% mtx_trns_q = mtx_trns./min(mtx_trns(:));
%plot(mtx_trns_q(:,3));

unique(mtx_trns(:,1)) ;   % 0.50746 transition exists here but not in column 2.
unique(mtx_trns(:,2));
% plot(unique(mtx_trns(:,1)),unique(mtx_trns(:,2)))

length(unique(mtx_trns(:,1)));
length(unique(mtx_trns(:,2)));

% NOTE: go back to stitching file to eliminate transitions that occur at
% the end of one file, where it turns to another.


% A = array of values from matrix transitions
% B = array of values within stdev from FRET state for each transition
% Lia = isamember(A,B) is a logical true/false statement
% this is going to be tricky because we need all the digists to match.

%A = [unique(mtx_trns(:,1))];
%B = [(0.37255-0.2):0.00001:(0.37255+0.2)];
%[C, iA, iB] = intersect(A,B);

roundFRET_i = round(mtx_trns(:,1),2);
roundFRET_f = round(mtx_trns(:,2),2);
roundFRET_states_i = unique(roundFRET_i);
roundFRET_states_f = unique(roundFRET_f);

%round_mtx_trns = [roundFRET_i,roundFRET_f,mtx_trns(:,3)]
round_mtx_trns = [roundFRET_i,roundFRET_f];

% Because of rounding errors caused by the meshgrid function, for some of
% the FRET states there are two values recorded by the unique function
% which we can consider the same state. So for each of these cases we will
% label them a and b. The error between these states is +/- 0.01, which is
% within the standard deviation of the HMM analysis done by vbFRET for each
% FRET state value. 

FRET_1 = roundFRET_states_i(2);
FRET_2_a = roundFRET_states_i(4);
FRET_2_b = roundFRET_states_i(5);
FRET_3_a = roundFRET_states_i(6);
FRET_3_b = roundFRET_states_i(7);
FRET_4_a = roundFRET_states_i(8);
FRET_4_b = roundFRET_states_i(9);
FRET_5 = roundFRET_states_i(10);

% define the possible FRET observations for each rate, considering both a
% and b values (defined above).
k12_FRET_a = [FRET_1, FRET_2_a];
k12_FRET_b = [FRET_1, FRET_2_b];

k21_FRET_a = [FRET_2_a, FRET_1];
k21_FRET_b = [FRET_2_b, FRET_1];

k13_FRET_a = [FRET_1, FRET_3_a];
k13_FRET_b = [FRET_1, FRET_3_b];

k31_FRET_a = [FRET_3_a, FRET_1];
k31_FRET_b = [FRET_3_b, FRET_1];

k14_FRET_a = [FRET_1, FRET_4_a];
k14_FRET_b = [FRET_1, FRET_4_b];

k41_FRET_a = [FRET_4_a, FRET_1];
k41_FRET_b = [FRET_4_b, FRET_1];

k15_FRET = [FRET_1, FRET_5];

k51_FRET = [FRET_5, FRET_1];

k23_FRET_aa = [FRET_2_a, FRET_3_a];
k23_FRET_ab  = [FRET_2_a, FRET_3_b];
k23_FRET_ba = [FRET_2_b, FRET_3_a];
k23_FRET_bb = [FRET_2_b, FRET_3_b];

k32_FRET_aa = [FRET_3_a, FRET_2_a];
k32_FRET_ab = [FRET_3_a, FRET_2_b];
k32_FRET_ba = [FRET_3_b, FRET_2_a];
k32_FRET_bb = [FRET_3_b, FRET_2_b];

k24_FRET_aa = [FRET_2_a, FRET_4_a];
k24_FRET_ab = [FRET_2_a, FRET_4_b];
k24_FRET_ba = [FRET_2_b, FRET_4_a];
k24_FRET_bb = [FRET_2_b, FRET_4_b];

k42_FRET_aa = [FRET_4_a, FRET_2_a];
k42_FRET_ab = [FRET_4_a, FRET_2_b];
k42_FRET_ba = [FRET_4_b, FRET_2_a];
k42_FRET_bb = [FRET_4_b, FRET_2_b];

k25_FRET_a = [FRET_2_a, FRET_5];
k25_FRET_b = [FRET_2_b, FRET_5];

k52_FRET_a = [FRET_5, FRET_2_a];
k52_FRET_b = [FRET_5, FRET_2_b];

k34_FRET_aa = [FRET_3_a, FRET_4_a];
k34_FRET_ab = [FRET_3_a, FRET_4_b];
k34_FRET_ba = [FRET_3_b, FRET_4_a];
k34_FRET_bb = [FRET_3_b, FRET_4_b];

k43_FRET_aa = [FRET_4_a, FRET_3_a];
k43_FRET_ab = [FRET_4_a, FRET_3_b];
k43_FRET_ba = [FRET_4_b, FRET_3_a];
k43_FRET_bb = [FRET_4_b, FRET_3_b];

k35_FRET_a = [FRET_3_a, FRET_5];
k35_FRET_b = [FRET_3_b, FRET_5];

k53_FRET_a = [FRET_5, FRET_3_a];
k53_FRET_b = [FRET_5, FRET_3_b];

k45_FRET_a = [FRET_4_a, FRET_5];
k45_FRET_b = [FRET_4_b, FRET_5];

k54_FRET_a = [FRET_5, FRET_4_a];
k54_FRET_b = [FRET_5, FRET_4_b];

%% Which rates exist in our data set?
% Define a column of the possible rates
rates_possible = [k12; k13; k14; k15; k21; k23; k24; k25; k31; k32; k34; k35; k41; k42; k43; k45; k51; k52; k53; k54];

% Ask the data set if the FRET transitions corresponding to each rate exist
%  in the set of rounded data.

% Does rate k12 exist? If so, write k12, if not k12 = 0.
if ismember(k12_FRET_a,round_mtx_trns,'rows')||ismember(k12_FRET_b,round_mtx_trns,'rows')
    k12;
    disp('k12 exists')
else
    k12 = 0;
    disp('k12 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k12]);
else
    rates = [k12];
end

% Does rate k13 exist? If so, write k13, if not k13 = 0.
if ismember(k13_FRET_a,round_mtx_trns,'rows')||ismember(k13_FRET_b,round_mtx_trns,'rows')
    k13;
    disp('k13 exists')
else
    k13 = 0;
    disp('k13 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k13]);
else
    rates = [k13];
end
            
% Does rate k14 exist? If so, write k14, if not k14 = 0.
if ismember(k14_FRET_a,round_mtx_trns,'rows')||ismember(k14_FRET_b,round_mtx_trns,'rows')
    k14;
    disp('k14 exists')
else
    k14 = 0;
    disp('k14 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k14]);
else
    rates = [k14];
end

% Does rate k15 exist? If so, write k15, if not k15 = 0.
if ismember(k15_FRET,round_mtx_trns,'rows')
    k15;
    disp('k15 exists')
else
    k15 = 0;
    disp('k15 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k15]);
else
    rates = [k15];
end

% Does rate k21 exist? If so, write k21, if not k21 = 0.
if ismember(k21_FRET_a,round_mtx_trns,'rows')||ismember(k21_FRET_b,round_mtx_trns,'rows')
    k21;
    disp('k21 exists')
else
    k21 = 0;
    disp('k21 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k21]);
else
    rates = [k21];
end

% Does rate k23 exist? If so, write k23, if not k23 = 0.
if ismember(k23_FRET_aa,round_mtx_trns,'rows')||ismember(k23_FRET_ab,round_mtx_trns,'rows')||ismember(k23_FRET_ba,round_mtx_trns,'rows')||ismember(k23_FRET_bb,round_mtx_trns,'rows')
    k23;
    disp('k23 exists')
else
    k23 = 0;
    disp('k23 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k23]);
else
    rates = [k23];
end

% Does rate k24 exist? If so, write k24, if not k24 = 0.
if ismember(k24_FRET_aa,round_mtx_trns,'rows')||ismember(k24_FRET_ab,round_mtx_trns,'rows')||ismember(k24_FRET_ba,round_mtx_trns,'rows')||ismember(k24_FRET_bb,round_mtx_trns,'rows')
    k24;
    disp('k24 exists')
else
    k24 = 0;
    disp('k24 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k24]);
else
    rates = [k24];
end

% Does rate k25 exist? If so, write k25, if not k25 = 0.
if ismember(k25_FRET_a,round_mtx_trns,'rows')||ismember(k25_FRET_b,round_mtx_trns,'rows')
    k25;
    disp('k25 exists')
else
    k25 = 0;
    disp('k25 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k25]);
else
    rates = [k25];
end

% Does rate k31 exist? If so, write k31, if not k31 = 0.
if ismember(k31_FRET_a,round_mtx_trns,'rows')||ismember(k31_FRET_b,round_mtx_trns,'rows')
    k31;
    disp('k31 exists')
else
    k31 = 0;
    disp('k31 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k31]);
else
    rates = [k31];
end

% Does rate k32 exist? If so, write k32, if not k32 = 0.
if ismember(k32_FRET_aa,round_mtx_trns,'rows')||ismember(k32_FRET_ab,round_mtx_trns,'rows')||ismember(k32_FRET_ba,round_mtx_trns,'rows')||ismember(k32_FRET_bb,round_mtx_trns,'rows')
    k32;
    disp('k32 exists')
else
    k32 = 0;
    disp('k32 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k32]);
else
    rates = [k32];
end


% Does rate k34 exist? If so, write k34, if not k34 = 0.
if ismember(k34_FRET_aa,round_mtx_trns,'rows')||ismember(k34_FRET_ab,round_mtx_trns,'rows')||ismember(k34_FRET_ba,round_mtx_trns,'rows')||ismember(k34_FRET_bb,round_mtx_trns,'rows')
    k34;
    disp('k34 exists')
else
    k34 = 0;
    disp('k34 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k34]);
else
    rates = [k34];
end

% Does rate k35 exist? If so, write k35, if not k35 = 0.
if ismember(k35_FRET_a,round_mtx_trns,'rows')||ismember(k35_FRET_b,round_mtx_trns,'rows')
    k35;
    disp('k35 exists')
else
    k35 = 0;
    disp('k35 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k35]);
else
    rates = [k35];
end

% Does rate k41 exist? If so, write k41, if not k41 = 0.
if ismember(k41_FRET_a,round_mtx_trns,'rows')||ismember(k41_FRET_b,round_mtx_trns,'rows')
    k41;
    disp('k41 exists')
else
    k41 = 0;
    disp('k41 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k41]);
else
    rates = [k41];
end

% Does rate k42 exist? If so, write k42, if not k42 = 0.
if ismember(k42_FRET_aa,round_mtx_trns,'rows')||ismember(k42_FRET_ab,round_mtx_trns,'rows')||ismember(k42_FRET_ba,round_mtx_trns,'rows')||ismember(k42_FRET_bb,round_mtx_trns,'rows')
    k42;
    disp('k42 exists')
else
    k42 = 0;
    disp('k42 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k42]);
else
    rates = [k42];
end

% Does rate k43 exist? If so, write k43, if not k43 = 0.
if ismember(k43_FRET_aa,round_mtx_trns,'rows')||ismember(k43_FRET_ab,round_mtx_trns,'rows')||ismember(k43_FRET_ba,round_mtx_trns,'rows')||ismember(k43_FRET_bb,round_mtx_trns,'rows')
    k43;
    disp('k43 exists')
else
    k43 = 0;
    disp('k43 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k43]);
else
    rates = [k43];
end

% Does rate k45 exist? If so, write k45, if not k45 = 0.
if ismember(k45_FRET_a,round_mtx_trns,'rows')||ismember(k45_FRET_b,round_mtx_trns,'rows')
    k45;
    disp('k45 exists')
else
    k45 = 0;
    disp('k45 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k45]);
else
    rates = [k45];
end

% Does rate k51 exist? If so, write k51, if not k51 = 0.
if ismember(k51_FRET,round_mtx_trns,'rows')
    k51;
    disp('k51 exists')
else
    k51 = 0;
    disp('k51 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k51]);
else
    rates = [k51];
end

% Does rate k52 exist? If so, write k52, if not k52 = 0.
if ismember(k52_FRET_a,round_mtx_trns,'rows')||ismember(k52_FRET_b,round_mtx_trns,'rows')
    k52;
    disp('k52 exists')
else
    k52 = 0;
    disp('k52 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k52]);
else
    rates = [k52];
end

% Does rate k53 exist? If so, write k53, if not k53 = 0.
if ismember(k53_FRET_a,round_mtx_trns,'rows')||ismember(k53_FRET_b,round_mtx_trns,'rows')
    k53;
    disp('k53 exists')
else
    k53 = 0;
    disp('k53 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k53]);
else
    rates = [k53];
end

% Does rate k54 exist? If so, write k54, if not k54 = 0.
if ismember(k54_FRET_a,round_mtx_trns,'rows')||ismember(k54_FRET_b,round_mtx_trns,'rows')
    k54;
    disp('k54 exists')
else
    k54 = 0;
    disp('k54 = 0')
end

if exist('rates')
    rates = vertcat(rates, [k54]);
else
    rates = [k54];
end

%%
% Create a matrix with a colum of possible rates, and the outcome of our
% above analysis about which rates are observed.
rates;
rates_obs = horzcat(rates_possible, rates)


save('dwellTDP_trim','fig1','rates_obs')


close all
clear all
clc
%%
% AUTHORS:  Claire Albrecht & Brett Israels
%
% NAME:     dwell2Prob.m
%
% PURPOSE:  Make histogram of dwell times and survival probability plot.
%
% INPUT:    dwell file (from pathDwell program) - only the stitched file
%
% OUTPUT:   Histogram of dwell times and plot of survival probability.
%           x-data: 
%           surv_prob_ij_norm, 
%
% NOTE 1:   the dwell times are in frames - not actual time
%           -> multiply by resolution to get time.
%
% CREATED:  July 2019
% MODIFIED: July 2019 - remove artificial transitions from stitching
%
%__________________________________________________________________________

%--------------------------------------------------------------------------
% Find dwell files - full dwell files
%--------------------------------------------------------------------------

<<<<<<< Updated upstream
%dwellFileName = '3p15mer_3.000000e-02trace_stitched_dwell.dat';
dwellFileName = '3p15mer_5.000000e-03trace_stitched_dwell.dat';

%res = 0.03;
%--------------------------------------------------------------------------
% Find resolution based on the filename (1-24-19)
%--------------------------------------------------------------------------
usres_str = dwellFileName(9:16);
%usres_str = datafiles(1).name(14:19);
usres = str2double(usres_str);
msres = usres*10^-3;%('i.e. 1000 usec = 1 msec')
res = usres*10^-6;%(i.e. 1 million useconds turn into 1 second)
disp(['Setting the resolution equal to ' num2str(usres) ' microseconds (' num2str(msres) ' msec) (' num2str(res) ' sec)']);

%--------------------------------------------------------------------------
% From dwell file

% dwelltimes = dlmread(dwellFileName);
% FRET_initial = dwelltimes(:,1);
% FRET_final = dwelltimes(:,2);
% dwell_time_frames = dwelltimes(:,3);
% 
% dwell_time = dwell_time_frames * res;
% 
% FRET_data = [FRET_initial, FRET_final, dwell_time];
% 
% FRET_states = unique(FRET_initial);
%--------------------------------------------------------------------------

% From chopped dwell file - variable
load('dwellInfoChopped.mat')
FRET_initial = dwellInfoChopped(:,1);
FRET_final = dwellInfoChopped(:,2);
dwell_time_frames = dwellInfoChopped(:,3);
=======
% dwellFileName = '3p15mer_3.000000e-02trace_stitched_dwell.dat';
% 
% res = 0.03;
% 
% dwelltimes = dlmread(dwellFileName);
% FRET_initial = dwelltimes(:,1);
% FRET_final = dwelltimes(:,2);
% dwell_time_frames = dwelltimes(:,3);
% fileNums = dwelltimes(:,4);

%%
%--------------------------------------------------------------------------
% Take data from the chopped file
%--------------------------------------------------------------------------

load('dwellInfoChopped.mat');
res = 0.03;
FRET_initial = dwellInfoChopped(:,1);
FRET_final = dwellInfoChopped(:,2);
dwell_time_frames = dwellInfoChopped(:,3);
fileNums = dwellInfoChopped(:,4);

%%
>>>>>>> Stashed changes

dwell_time = dwell_time_frames * res;

FRET_data = [FRET_initial, FRET_final, dwell_time, fileNums];

FRET_states = nonzeros(unique(FRET_initial));


%% 
%--------------------------------------------------------------------------
% Eliminate end of file transitions:
%--------------------------------------------------------------------------

fileEndLoc = find(FRET_data(:,4)==0);

% At file ends:
FRET_initial(fileEndLoc,1) = 0;    %   Set initial FRET to zero
FRET_final(fileEndLoc,2) = 0;      %   Set final FRET to zero

FRET_data(fileEndLoc,[1, 2]) = 0;  % Write this into the FRET_data variable
                                   % so they do not register as transitions for
                                   % the rest of the calculation.

<<<<<<< Updated upstream
FRET_states = nonzeros(unique(FRET_initial));
=======

>>>>>>> Stashed changes

%%
%--------------------------------------------------------------------------
% For transition 1 to 2:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(2)) > 0
    disp("k12 exists")
    rate_ID = 'k12';
    trns_12 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(2), :);
    dwellTimes_12 = trns_12(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_12);
    NBins = max(dwellTimes_12)/res + 1; 
    % Use these to define the edges of the histogram bins
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 1 --> 2")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    %counts_12 = histcounts(dwellTimes_12,NBins);
    counts_12 = histcounts(dwellTimes_12,timeBarEdges);
    timePlot_12 = timeBarEdges(1:end-1);
    bar(timePlot_12,counts_12)
    
    N_tot = sum(counts_12);
    temp_tot = 0;
    surv_prob_12 = [];
    
    i = 1;
    while i <= numel(counts_12)
        
        temp = counts_12(i);
        temp_tot = temp_tot + temp;
        
        if counts_12(i) == 0    % If this is satisfied, we do not want to collect the point.
            i = i + 1;          % Continue does not automatically increase i, so need this or else stuck in infinite loop.
            continue            % This tells the program to continue to next step of loop.
        end                     % This section of code eliminates straight lines between points - need scatter plot.
        
        surv_prob_12 = vertcat(surv_prob_12, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
%     time = [];
%     for i = 1:numel(counts_12)
%     time = [time, timeBarEdges(counts_12(i) > 0)];       % This will select the times for which the condition is met.
%     end
    %time = timeBarEdges(counts_12 > 0);
    
    time_counts_12 = find(counts_12 > 0);
    time = timeBarEdges(time_counts_12);
    prob = surv_prob_12 / max(surv_prob_12);    % Normalize the probability.
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 1 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);              % This is just for saving convenience, want columns not rows.
    
    save('survProb_12.mat','time','prob','rate_ID')
    
else
    disp("k12 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 1 to 3:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(3)) > 0
    disp("k13 exists")
    rate_ID = 'k13';
    trns_13 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(3), :);
    dwellTimes_13 = trns_13(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_13);
    NBins = max(dwellTimes_13)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition.
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 1 --> 3")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_13 = histcounts(dwellTimes_13,timeBarEdges);
    timePlot_13 = timeBarEdges(1:end-1);
    bar(timePlot_13,counts_13)
    
    N_tot = sum(counts_13);
    temp_tot = 0;
    surv_prob_13 = [];
    
    i = 1;
    while i <= numel(counts_13)
        
        temp = counts_13(i);
        temp_tot = temp_tot + temp;
        
        if counts_13(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_13 = vertcat(surv_prob_13, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
    %time = timeBarEdges( counts_13 > 0 );
    time_counts_13 = find(counts_13 > 0);
    time = timeBarEdges(time_counts_13);
    prob = surv_prob_13 / max(surv_prob_13);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 1 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_13.mat','time','prob','rate_ID')
    
else
    disp("k13 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 1 to 4: 
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(4)) > 0
    disp("k14 exists")
    rate_ID = 'k14';
    trns_14 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(4), :);
    dwellTimes_14 = trns_14(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_14);
    NBins = max(dwellTimes_14)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition.
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 1 --> 4")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_14 = histcounts(dwellTimes_14,timeBarEdges);
    timePlot_14 = timeBarEdges(1:end-1);
    bar(timePlot_14,counts_14)
    
    N_tot = sum(counts_14);
    temp_tot = 0;
    surv_prob_14 = [];
    
    i = 1;
    while i <= numel(counts_14)
        
        temp = counts_14(i);
        temp_tot = temp_tot + temp;
        
        if counts_14(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_14 = vertcat(surv_prob_14, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  %  time = timeBarEdges( counts_14 > 0 );
    time_counts_14 = find(counts_14 > 0);
    time = timeBarEdges(time_counts_14);  
    prob = surv_prob_14 / max(surv_prob_14);

    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 1 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
  
    time = transpose(time);
    
    save('survProb_14.mat','time','prob','rate_ID')
    
else
    disp("k14 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 1 to 5:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(5)) > 0
    disp("k15 exists")
    rate_ID = 'k15';
    trns_15 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_15 = trns_15(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_15);
    NBins = max(dwellTimes_15)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 1 --> 5")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_15 = histcounts(dwellTimes_15,timeBarEdges);
    timePlot_15 = timeBarEdges(1:end-1);
    bar(timePlot_15,counts_15)
    
    N_tot = sum(counts_15);
    temp_tot = 0;
    surv_prob_15 = [];
    
    i = 1;
    while i <= numel(counts_15)
        
        temp = counts_15(i);
        temp_tot = temp_tot + temp;
        
        if counts_15(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_15 = vertcat(surv_prob_15, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
%    time = timeBarEdges( counts_15 > 0 );
    time_counts_15 = find(counts_15 > 0);
    time = timeBarEdges(time_counts_15);
    prob = surv_prob_15 / max(surv_prob_15);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 1 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_15.mat','time','prob','rate_ID')
    
else
    disp("k15 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 1:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(1) ) > 0
    disp("k21 exists")
    rate_ID = 'k21';
    trns_21 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_21 = trns_21(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_21);
    NBins = max(dwellTimes_21)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 2 --> 1")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_21 = histcounts(dwellTimes_21,timeBarEdges);
    timePlot_21 = timeBarEdges(1:end-1);
    bar(timePlot_21,counts_21)
    
    N_tot = sum(counts_21);
    temp_tot = 0;
    surv_prob_21 = [];
    
    i = 1;
    while i <= numel(counts_21)
        
        temp = counts_21(i);
        temp_tot = temp_tot + temp;
        
        if counts_21(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_21 = vertcat(surv_prob_21, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  %  time = timeBarEdges( counts_21 > 0 );
    time_counts_21 = find(counts_21 > 0);
    time = timeBarEdges(time_counts_21);
    prob = surv_prob_21 / max(surv_prob_21);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 2 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_21.mat','time','prob','rate_ID')
    
else
    disp("k21 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 3:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(3) ) > 0
    disp("k23 exists")
    rate_ID = 'k23';
    trns_23 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(3), :);
    dwellTimes_23 = trns_23(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_23);
    NBins = max(dwellTimes_23)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 2 --> 3")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_23 = histcounts(dwellTimes_23,timeBarEdges);
    timePlot_23 = timeBarEdges(1:end-1);
    bar(timePlot_23,counts_23)
    
    N_tot = sum(counts_23);
    temp_tot = 0;
    surv_prob_23 = [];
    
    i = 1;
    while i <= numel(counts_23)
        
        temp = counts_23(i);
        temp_tot = temp_tot + temp;
        
        if counts_23(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_23 = vertcat(surv_prob_23, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  %  time = timeBarEdges( counts_23 > 0 );
    time_counts_23 = find(counts_23 > 0);
    time = timeBarEdges(time_counts_23);
    prob = surv_prob_23 / max(surv_prob_23);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 2 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_23.mat','time','prob','rate_ID')
    
else
    disp("k23 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 4:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(4) ) > 0
    disp("k24 exists")
    rate_ID = 'k24';
    trns_24 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(4), :);
    dwellTimes_24 = trns_24(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_24);
    NBins = max(dwellTimes_24)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 2 --> 4")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_24 = histcounts(dwellTimes_24,timeBarEdges);
    timePlot_24 = timeBarEdges(1:end-1);
    bar(timePlot_24,counts_24)
    
    N_tot = sum(counts_24);
    temp_tot = 0;
    surv_prob_24 = [];
    
    i = 1;
    while i <= numel(counts_24)
        
        temp = counts_24(i);
        temp_tot = temp_tot + temp;
        
        if counts_24(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_24 = vertcat(surv_prob_24, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  %  time = timeBarEdges( counts_24 > 0 );
    time_counts_24 = find(counts_24 > 0);
    time = timeBarEdges(time_counts_24);
    prob = surv_prob_24 / max(surv_prob_24);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 2 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_24.mat','time','prob','rate_ID')
    
else
    disp("k24 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 5:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(5) ) > 0
    disp("k25 exists")
    rate_ID = 'k25';
    trns_25 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_25 = trns_25(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_25);
    NBins = max(dwellTimes_25)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition.
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 2 --> 5")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_25 = histcounts(dwellTimes_25,timeBarEdges);
    timePlot_25 = timeBarEdges(1:end-1);
    bar(timePlot_25,counts_25)
    
    N_tot = sum(counts_25);
    temp_tot = 0;
    surv_prob_25 = [];
    
    i = 1;
    while i <= numel(counts_25)
        
        temp = counts_25(i);
        temp_tot = temp_tot + temp;
        
        if counts_25(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_25 = vertcat(surv_prob_25, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
 %   time = timeBarEdges( counts_25 > 0 );
    time_counts_25 = find(counts_25 > 0);
    time = timeBarEdges(time_counts_25);
    prob = surv_prob_25 / max(surv_prob_25);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 2 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_25.mat','time','prob', 'rate_ID')
    
else
    disp("k25 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 1:
%--------------------------------------------------------------------------
if sum ( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(1) ) > 0
    disp("k31 exists")
    rate_ID = 'k31';
    trns_31 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_31 = trns_31(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_31);
    NBins = max(dwellTimes_31)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 3 --> 1")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_31 = histcounts(dwellTimes_31,timeBarEdges);
    timePlot_31 = timeBarEdges(1:end-1);
    bar(timePlot_31,counts_31)
    
    N_tot = sum(counts_31);
    temp_tot = 0;
    surv_prob_31 = [];
    
    i = 1;
    while i <= numel(counts_31)
        
        temp = counts_31(i);
        temp_tot = temp_tot + temp;
        
        if counts_31(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_31 = vertcat(surv_prob_31, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
%    time = timeBarEdges( counts_31 > 0 );
    time_counts_31 = find(counts_31 > 0);
    time = timeBarEdges(time_counts_31);
    prob = surv_prob_31 / max(surv_prob_31);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 3 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_31.mat','time','prob','rate_ID')
    
else
    disp("k31 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 2:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(2) ) > 0
    disp("k32 exists")
    rate_ID = 'k32';
    trns_32 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(2), :);
    dwellTimes_32 = trns_32(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_32);
    NBins = max(dwellTimes_32)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 3 --> 2")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_32 = histcounts(dwellTimes_32,timeBarEdges);
    timePlot_32 = timeBarEdges(1:end-1);
    bar(timePlot_32,counts_32)
    
    N_tot = sum(counts_32);
    temp_tot = 0;
    surv_prob_32 = [];
    
    i = 1;
    while i <= numel(counts_32)
        
        temp = counts_32(i);
        temp_tot = temp_tot + temp;
        
        if counts_32(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_32 = vertcat(surv_prob_32, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
%    time = timeBarEdges( counts_32 > 0 );
    time_counts_32 = find(counts_32 > 0);
    time = timeBarEdges(time_counts_32);
    prob = surv_prob_32 / max(surv_prob_32);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 3 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_32.mat','time','prob','rate_ID')
    
else
    disp("k32 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 4:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(4) ) > 0
    disp("k34 exists")
    rate_ID = 'k34';
    trns_34 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(4), :);
    dwellTimes_34 = trns_34(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_34);
    NBins = max(dwellTimes_34)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 3 --> 4")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_34 = histcounts(dwellTimes_34,timeBarEdges);
    timePlot_34 = timeBarEdges(1:end-1);
    bar(timePlot_34,counts_34)
    
    N_tot = sum(counts_34);
    temp_tot = 0;
    surv_prob_34 = [];
    
    i = 1;
    while i <= numel(counts_34)
        
        temp = counts_34(i);
        temp_tot = temp_tot + temp;
        
        if counts_34(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_34 = vertcat(surv_prob_34, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  %  time = timeBarEdges( counts_34 > 0 );
    time_counts_34 = find(counts_34 > 0);
    time = timeBarEdges(time_counts_34);
    prob = surv_prob_34 / max(surv_prob_34);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 3 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_34.mat','time','prob','rate_ID')
    
else
    disp("k34 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 5:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(5) ) > 0
    disp("k35 exists")
    rate_ID = 'k35';
    trns_35 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_35 = trns_35(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_35);
    NBins = max(dwellTimes_35)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 3 --> 5")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_35 = histcounts(dwellTimes_35,timeBarEdges);
    timePlot_35 = timeBarEdges(1:end-1);
    bar(timePlot_35,counts_35)
    
    N_tot = sum(counts_35);
    temp_tot = 0;
    surv_prob_35 = [];
    
    i = 1;
    while i <= numel(counts_35)
        
        temp = counts_35(i);
        temp_tot = temp_tot + temp;
        
        if counts_35(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_35 = vertcat(surv_prob_35, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
   % time = timeBarEdges( counts_35 > 0 );
    time_counts_35 = find(counts_35 > 0);
    time = timeBarEdges(time_counts_35);
    prob = surv_prob_35 / max(surv_prob_35);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 3 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_35.mat','time','prob','rate_ID')
    
else
    disp("k35 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 1:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(1) ) > 0
    disp("k41 exists")
    rate_ID = 'k41';
    trns_41 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_41 = trns_41(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_41);
    NBins = max(dwellTimes_41)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 4 --> 1")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_41 = histcounts(dwellTimes_41,timeBarEdges);
    timePlot_41 = timeBarEdges(1:end-1);
    bar(timePlot_41,counts_41)
    
    N_tot = sum(counts_41);
    temp_tot = 0;
    surv_prob_41 = [];
    
    i = 1;
    while i <= numel(counts_41)
        
        temp = counts_41(i);
        temp_tot = temp_tot + temp;
        
        if counts_41(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_41 = vertcat(surv_prob_41, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  %  time = timeBarEdges( counts_41 > 0 );
    time_counts_41 = find(counts_41 > 0);
    time = timeBarEdges(time_counts_41);
    prob = surv_prob_41 / max(surv_prob_41);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 4 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_41.mat','time','prob','rate_ID')
    
else
    disp("k41 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 2:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(2) ) > 0
    disp("k42 exists")
    rate_ID = 'k42';
    trns_42 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(2), :);
    dwellTimes_42 = trns_42(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_42);
    NBins = max(dwellTimes_42)/res + 1;  
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 4 --> 2")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_42 = histcounts(dwellTimes_42,timeBarEdges);
    timePlot_42 = timeBarEdges(1:end-1);
    bar(timePlot_42,counts_42)
    
    N_tot = sum(counts_42);
    temp_tot = 0;
    surv_prob_42 = [];
    
    i = 1;
    while i <= numel(counts_42)
        
        temp = counts_42(i);
        temp_tot = temp_tot + temp;
        
        if counts_42(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_42 = vertcat(surv_prob_42, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
 %   time = timeBarEdges( counts_42 > 0 );
    time_counts_42 = find(counts_42 > 0);
    time = timeBarEdges(time_counts_42);
    prob = surv_prob_42 / max(surv_prob_42);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 4 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_42.mat','time','prob','rate_ID')
    
else
    disp("k42 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 3:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(3)) > 0
    disp("k43 exists")
    rate_ID = 'k43';
    trns_43 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(3), :);
    dwellTimes_43 = trns_43(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_43);
    NBins = max(dwellTimes_43)/res + 1;  
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 4 --> 3")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_43 = histcounts(dwellTimes_43,timeBarEdges);
    timePlot_43 = timeBarEdges(1:end-1);
    bar(timePlot_43,counts_43)
    
    N_tot = sum(counts_43);
    temp_tot = 0;
    surv_prob_43 = [];
    
    i = 1;
    while i <= numel(counts_43)
        
        temp = counts_43(i);
        temp_tot = temp_tot + temp;
        
        if counts_43(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_43 = vertcat(surv_prob_43, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  %  time = timeBarEdges( counts_43 > 0 );
    time_counts_43 = find(counts_43 > 0);
    time = timeBarEdges(time_counts_43);
    prob = surv_prob_43 / max(surv_prob_43);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 4 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_43.mat','time','prob','rate_ID')
    
else
    disp("k43 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 5:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(5)) > 0
    disp("k45 exists")
    rate_ID = 'k45';
    trns_45 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_45 = trns_45(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_45);
    NBins = max(dwellTimes_45)/res + 1;
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 4 --> 5")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_45 = histcounts(dwellTimes_45,timeBarEdges);
    timePlot_45 = timeBarEdges(1:end-1);
    bar(timePlot_45,counts_45)
    
    N_tot = sum(counts_45);
    temp_tot = 0;
    surv_prob_45 = [];
    
    i = 1;
    while i <= numel(counts_45)
        
        temp = counts_45(i);
        temp_tot = temp_tot + temp;
        
        if counts_45(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_45 = vertcat(surv_prob_45, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
   % time = timeBarEdges( counts_45 > 0 );
    time_counts_45 = find(counts_45 > 0);
    time = timeBarEdges(time_counts_45);
    prob = surv_prob_45 / max(surv_prob_45);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 4 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_45.mat','time','prob','rate_ID')
    
else
    disp("k45 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 1:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(1)) > 0
    disp("k51 exists")
    rate_ID = 'k51';
    trns_51 = FRET_data( FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_51 = trns_51(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_51);
    NBins = max(dwellTimes_51)/res + 1; 
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 5 --> 1")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_51 = histcounts(dwellTimes_51,timeBarEdges);
    timePlot_51 = timeBarEdges(1:end-1);
    bar(timePlot_51,counts_51)
    
    N_tot = sum(counts_51);
    temp_tot = 0;
    surv_prob_51 = [];
    
    i = 1;
    while i <= numel(counts_51)
        
        temp = counts_51(i);
        temp_tot = temp_tot + temp;
        
        if counts_51(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_51 = vertcat(surv_prob_51, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
 %   time = timeBarEdges( counts_51 > 0 );
    time_counts_51 = find(counts_51 > 0);
    time = timeBarEdges(time_counts_51);
    prob = surv_prob_51 / max(surv_prob_51);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 5 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_51.mat','time','prob','rate_ID')
    
else
    disp("k51 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 2:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(2)) > 0
    disp("k52 exists")
    rate_ID = 'k52';
    trns_52 = FRET_data( FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(2), :);
    dwellTimes_52 = trns_52(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_52);
    NBins = max(dwellTimes_52)/res + 1; 
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 5 --> 2")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_52 = histcounts(dwellTimes_52,timeBarEdges);
    timePlot_52 = timeBarEdges(1:end-1);
    bar(timePlot_52,counts_52)
    
    N_tot = sum(counts_52);
    temp_tot = 0;
    surv_prob_52 = [];
    
    i = 1;
    while i <= numel(counts_52)
        
        temp = counts_52(i);
        temp_tot = temp_tot + temp;
        
        if counts_52(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_52 = vertcat(surv_prob_52, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
   % time = timeBarEdges( counts_52 > 0 );
    time_counts_52 = find(counts_52 > 0);
    time = timeBarEdges(time_counts_52);
    prob = surv_prob_52 / max(surv_prob_52);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 5 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_52.mat','time','prob','rate_ID')
    
else
    disp("k52 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 3:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(3)) > 0
    disp("k53 exists")
    rate_ID = 'k53';
    trns_53 = FRET_data( FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(3), :);
    dwellTimes_53 = trns_53(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_53);
    NBins = max(dwellTimes_53)/res + 1; 
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 5 --> 3")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_53 = histcounts(dwellTimes_53,timeBarEdges);
    timePlot_53 = timeBarEdges(1:end-1);
    bar(timePlot_53,counts_53)
    
    N_tot = sum(counts_53);
    temp_tot = 0;
    surv_prob_53 = [];
    
    i = 1;
    while i <= numel(counts_53)
        
        temp = counts_53(i);
        temp_tot = temp_tot + temp;
        
        if counts_53(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_53 = vertcat(surv_prob_53, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
   % time = timeBarEdges( counts_53 > 0 );
    time_counts_53 = find(counts_53 > 0);
    time = timeBarEdges(time_counts_53);
    prob = surv_prob_53 / max(surv_prob_53);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 5 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_53.mat','time','prob','rate_ID')
    
else
    disp("k53 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 4:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(4)) > 0
    disp("k54 exists")
    rate_ID = 'k54';
    trns_54 = FRET_data( FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(4), :);
    dwellTimes_54 = trns_54(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_54);
    NBins = max(dwellTimes_54)/res + 1; 
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition
    figure
    subplot(1,2,1)
    title("Histogram of dwell times: 5 --> 4")
    xlabel("Dwell times (s)")
    ylabel("Counts")
    hold on
   
    counts_54 = histcounts(dwellTimes_54,timeBarEdges);
    timePlot_54 = timeBarEdges(1:end-1);
    bar(timePlot_54,counts_54);
    
    N_tot = sum(counts_54);
    temp_tot = 0;
    surv_prob_54 = [];
    
    i = 1;
    while i <= numel(counts_54)
        
        temp = counts_54(i);
        temp_tot = temp_tot + temp;
        
        if counts_54(i) == 0 
            i = i + 1;
            continue
        end
        
        surv_prob_54 = vertcat(surv_prob_54, (1 - (temp_tot / N_tot)));
        i = i + 1;
    end
    
  % time = timeBarEdges( counts_54 > 0 );
    time_counts_54 = find(counts_54 > 0);
    time = timeBarEdges(time_counts_54);
    prob = surv_prob_54 / max(surv_prob_54);
    
    subplot(1,2,2)
    scatter(time, prob);
    title("Survival probability: 5 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+res])
    hold off
    
    time = transpose(time);
    
    save('survProb_54.mat','time','prob','rate_ID')
    
else
    disp("k54 is zero")
end

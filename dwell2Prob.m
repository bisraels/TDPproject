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
% MODIFIED:
%
% -------------------------------------------------------------------------

% Find dwell files

dwellFileName = '3p15mer_3.000000e-02trace_stitched_dwell.dat';

res = 0.03;

dwelltimes = dlmread(dwellFileName);
FRET_initial = dwelltimes(:,1);
FRET_final = dwelltimes(:,2);
dwell_time_frames = dwelltimes(:,3);

dwell_time = dwell_time_frames * res;

FRET_data = [FRET_initial, FRET_final, dwell_time];

FRET_states = unique(FRET_initial);

% Define the FRET state transitions as they will appear in the data set.
FRET_12 = [FRET_states(1), FRET_states(2)];
FRET_13 = [FRET_states(1), FRET_states(3)];
FRET_14 = [FRET_states(1), FRET_states(4)];
FRET_15 = [FRET_states(1), FRET_states(5)];
FRET_21 = [FRET_states(2), FRET_states(1)];
FRET_23 = [FRET_states(2), FRET_states(3)];
FRET_24 = [FRET_states(2), FRET_states(4)];
FRET_25 = [FRET_states(2), FRET_states(5)];
FRET_31 = [FRET_states(3), FRET_states(1)];
FRET_32 = [FRET_states(3), FRET_states(2)];
FRET_34 = [FRET_states(3), FRET_states(4)];
FRET_35 = [FRET_states(3), FRET_states(5)];
FRET_41 = [FRET_states(4), FRET_states(1)];
FRET_42 = [FRET_states(4), FRET_states(2)];
FRET_43 = [FRET_states(4), FRET_states(3)];
FRET_45 = [FRET_states(4), FRET_states(5)];
FRET_51 = [FRET_states(5), FRET_states(1)];
FRET_52 = [FRET_states(5), FRET_states(2)];
FRET_53 = [FRET_states(5), FRET_states(3)];
FRET_54 = [FRET_states(5), FRET_states(4)];

%%
%--------------------------------------------------------------------------
% For transition 1 to 2:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(2)) > 0
    disp("k12 exists")
    trns_12 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(2), :);
    dwellTimes_12 = trns_12(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_12);
    NBins = max(dwellTimes_12)/res + 1;  % Not sure if this is right
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
    
    for i = 1:numel(counts_12)
        
        temp = counts_12(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_12 = vertcat(surv_prob_12, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_12_norm = surv_prob_12 / max(surv_prob_12);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_12, surv_prob_12_norm);
    title("Survival probability: 1 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_12 = transpose(timePlot_12);
    
    save('survProb_12.mat','timeVector_12','surv_prob_12_norm')
    
else
    disp("k12 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 1 to 3:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(3)) > 0
    disp("k13 exists")
    trns_13 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(3), :);
    dwellTimes_13 = trns_13(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_13);
    NBins = max(dwellTimes_13)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_13)
        
        temp = counts_13(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_13 = vertcat(surv_prob_13, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_13_norm = surv_prob_13 / max(surv_prob_13);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_13, surv_prob_13_norm);
    title("Survival probability: 1 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_13 = transpose(timePlot_13);
    
    save('survProb_13.mat','timeVector_13','surv_prob_13_norm')
    
else
    disp("k13 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 1 to 4:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(4)) > 0
    disp("k14 exists")
    trns_14 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(4), :);
    dwellTimes_14 = trns_14(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_14);
    NBins = max(dwellTimes_14)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_14)
        
        temp = counts_14(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_14 = vertcat(surv_prob_14, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_14_norm = surv_prob_14 / max(surv_prob_14);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_14, surv_prob_14_norm);
    title("Survival probability: 1 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_14 = transpose(timePlot_14);
    
    save('survProb_14.mat','timeVector_14','surv_prob_14_norm')
    
else
    disp("k14 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 1 to 5:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(5)) > 0
    disp("k15 exists")
    trns_15 = FRET_data( FRET_data(:,1) == FRET_states(1) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_15 = trns_15(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_15);
    NBins = max(dwellTimes_15)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_15)
        
        temp = counts_15(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_15 = vertcat(surv_prob_15, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_15_norm = surv_prob_15 / max(surv_prob_15);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_15, surv_prob_15_norm);
    title("Survival probability: 1 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_15 = transpose(timePlot_15);
    
    save('survProb_15.mat','timeVector_15','surv_prob_15_norm')
    
else
    disp("k15 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 1:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(1) ) > 0
    disp("k21 exists")
    trns_21 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_21 = trns_21(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_21);
    NBins = max(dwellTimes_21)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_21)
        
        temp = counts_21(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_21 = vertcat(surv_prob_21, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_21_norm = surv_prob_21 / max(surv_prob_21);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_21, surv_prob_21_norm);
    title("Survival probability: 2 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_21 = transpose(timePlot_21);
    
    save('survProb_21.mat','timeVector_21','surv_prob_21_norm')
    
else
    disp("k21 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 3:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(3) ) > 0
    disp("k23 exists")
    trns_23 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(3), :);
    dwellTimes_23 = trns_23(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_23);
    NBins = max(dwellTimes_23)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_23)
        
        temp = counts_23(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_23 = vertcat(surv_prob_23, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_23_norm = surv_prob_23 / max(surv_prob_23);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_23, surv_prob_23_norm);
    title("Survival probability: 2 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_23 = transpose(timePlot_23);
    
    save('survProb_23.mat','timeVector_23','surv_prob_23_norm')
    
else
    disp("k23 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 4:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(4) ) > 0
    disp("k24 exists")
    trns_24 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(4), :);
    dwellTimes_24 = trns_24(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_24);
    NBins = max(dwellTimes_24)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_24)
        
        temp = counts_24(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_24 = vertcat(surv_prob_24, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_24_norm = surv_prob_24 / max(surv_prob_24);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_24, surv_prob_24_norm);
    title("Survival probability: 2 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_24 = transpose(timePlot_24);
    
    save('survProb_24.mat','timeVector_24','surv_prob_24_norm')
    
else
    disp("k24 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 2 to 5:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(5) ) > 0
    disp("k25 exists")
    trns_25 = FRET_data( FRET_data(:,1) == FRET_states(2) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_25 = trns_25(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_25);
    NBins = max(dwellTimes_25)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_25)
        
        temp = counts_25(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_25 = vertcat(surv_prob_25, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_25_norm = surv_prob_25 / max(surv_prob_25);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_25, surv_prob_25_norm);
    title("Survival probability: 2 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_25 = transpose(timePlot_25);
    
    save('survProb_25.mat','timeVector_25','surv_prob_25_norm')
    
else
    disp("k25 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 1:
%--------------------------------------------------------------------------
if sum ( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(1) ) > 0
    disp("k31 exists")
    trns_31 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_31 = trns_31(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_31);
    NBins = max(dwellTimes_31)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 1 --> 2.
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
    
    for i = 1:numel(counts_31)
        
        temp = counts_31(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_31 = vertcat(surv_prob_31, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_31_norm = surv_prob_31 / max(surv_prob_31);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_31, surv_prob_31_norm);
    title("Survival probability: 3 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_31 = transpose(timePlot_31);
    
    save('survProb_31.mat','timeVector_31','surv_prob_31_norm')
    
else
    disp("k31 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 2:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(2) ) > 0
    disp("k32 exists")
    trns_32 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(2), :);
    dwellTimes_32 = trns_32(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_32);
    NBins = max(dwellTimes_32)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 3 --> 2.
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
    
    for i = 1:numel(counts_32)
        
        temp = counts_32(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_32 = vertcat(surv_prob_32, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_32_norm = surv_prob_32 / max(surv_prob_32);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_32, surv_prob_32_norm);
    title("Survival probability: 3 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_32 = transpose(timePlot_32);
    
    save('survProb_32.mat','timeVector_32','surv_prob_32_norm')
    
else
    disp("k32 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 4:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(4) ) > 0
    disp("k34 exists")
    trns_34 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(4), :);
    dwellTimes_34 = trns_34(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_34);
    NBins = max(dwellTimes_34)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 3 --> 2.
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
    
    for i = 1:numel(counts_34)
        
        temp = counts_34(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_34 = vertcat(surv_prob_34, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_34_norm = surv_prob_34 / max(surv_prob_34);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_34, surv_prob_34_norm);
    title("Survival probability: 3 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_34 = transpose(timePlot_34);
    
    save('survProb_34.mat','timeVector_34','surv_prob_34_norm')
    
else
    disp("k34 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 3 to 5:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(5) ) > 0
    disp("k35 exists")
    trns_35 = FRET_data( FRET_data(:,1) == FRET_states(3) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_35 = trns_35(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_35);
    NBins = max(dwellTimes_35)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 3 --> 5.
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
    
    for i = 1:numel(counts_35)
        
        temp = counts_35(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_35 = vertcat(surv_prob_35, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_35_norm = surv_prob_35 / max(surv_prob_35);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_35, surv_prob_35_norm);
    title("Survival probability: 3 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_35 = transpose(timePlot_35);
    
    save('survProb_35.mat','timeVector_35','surv_prob_35_norm')
    
else
    disp("k35 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 1:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(1) ) > 0
    disp("k41 exists")
    trns_41 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_41 = trns_41(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_41);
    NBins = max(dwellTimes_41)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 4 --> 1.
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
    
    for i = 1:numel(counts_41)
        
        temp = counts_41(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_41 = vertcat(surv_prob_41, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_41_norm = surv_prob_41 / max(surv_prob_41);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_41, surv_prob_41_norm);
    title("Survival probability: 4 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_41 = transpose(timePlot_41);
    
    save('survProb_41.mat','timeVector_41','surv_prob_41_norm')
    
else
    disp("k41 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 2:
%--------------------------------------------------------------------------
if sum( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(2) ) > 0
    disp("k42 exists")
    trns_42 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(2), :);
    dwellTimes_42 = trns_42(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_42);
    NBins = max(dwellTimes_42)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 4 --> 1.
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
    
    for i = 1:numel(counts_42)
        
        temp = counts_42(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_42 = vertcat(surv_prob_42, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_42_norm = surv_prob_42 / max(surv_prob_42);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_42, surv_prob_42_norm);
    title("Survival probability: 4 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_42 = transpose(timePlot_42);
    
    save('survProb_42.mat','timeVector_42','surv_prob_42_norm')
    
else
    disp("k42 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 3:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(3)) > 0
    disp("k43 exists")
    trns_43 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(3), :);
    dwellTimes_43 = trns_43(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_43);
    NBins = max(dwellTimes_43)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 4 --> 3.
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
    
    for i = 1:numel(counts_43)
        
        temp = counts_43(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_43 = vertcat(surv_prob_43, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_43_norm = surv_prob_43 / max(surv_prob_43);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_43, surv_prob_43_norm);
    title("Survival probability: 4 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_43 = transpose(timePlot_43);
    
    save('survProb_43.mat','timeVector_43','surv_prob_43_norm')
    
else
    disp("k43 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 4 to 5:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(5)) > 0
    disp("k45 exists")
    trns_45 = FRET_data( FRET_data(:,1) == FRET_states(4) & FRET_data(:,2) == FRET_states(5), :);
    dwellTimes_45 = trns_45(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_45);
    NBins = max(dwellTimes_45)/res + 1;  % Not sure if this is right
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 4 --> 3.
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
    
    for i = 1:numel(counts_45)
        
        temp = counts_45(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_45 = vertcat(surv_prob_45, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_45_norm = surv_prob_45 / max(surv_prob_45);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_45, surv_prob_45_norm);
    title("Survival probability: 4 --> 5")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_45 = transpose(timePlot_45);
    
    save('survProb_45.mat','timeVector_45','surv_prob_45_norm')
    
else
    disp("k45 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 1:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(1)) > 0
    disp("k51 exists")
    trns_51 = FRET_data( FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(1), :);
    dwellTimes_51 = trns_51(:,3);
    
    % Define lower and upper  bound for the FRET data.
    time_LB = 0;
    time_UB = max(dwellTimes_51);
    NBins = max(dwellTimes_51)/res + 1; 
    timeBarEdges  = [time_LB:(time_UB - time_LB)/(NBins -1):time_UB];
    
    % Make histogram of dwell times corresponding to transition 4 --> 3.
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
    
    for i = 1:numel(counts_51)
        
        temp = counts_51(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_51 = vertcat(surv_prob_51, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_51_norm = surv_prob_51 / max(surv_prob_51);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_51, surv_prob_51_norm);
    title("Survival probability: 5 --> 1")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_51 = transpose(timePlot_51);
    
    save('survProb_51.mat','timeVector_51','surv_prob_51_norm')
    
else
    disp("k51 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 2:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(2)) > 0
    disp("k52 exists")
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
    
    for i = 1:numel(counts_52)
        
        temp = counts_52(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_52 = vertcat(surv_prob_52, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_52_norm = surv_prob_52 / max(surv_prob_52);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_52, surv_prob_52_norm);
    title("Survival probability: 5 --> 2")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_52 = transpose(timePlot_52);
    
    save('survProb_52.mat','timeVector_52','surv_prob_52_norm')
    
else
    disp("k52 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 3:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(3)) > 0
    disp("k53 exists")
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
    
    for i = 1:numel(counts_53)
        
        temp = counts_53(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_53 = vertcat(surv_prob_53, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_53_norm = surv_prob_53 / max(surv_prob_53);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_53, surv_prob_53_norm);
    title("Survival probability: 5 --> 3")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_53 = transpose(timePlot_53);
    
    save('survProb_53.mat','timeVector_53','surv_prob_53_norm')
    
else
    disp("k53 is zero")
end

%%
%--------------------------------------------------------------------------
% For transition 5 to 4:
%--------------------------------------------------------------------------
if sum(FRET_data(:,1) == FRET_states(5) & FRET_data(:,2) == FRET_states(4)) > 0
    disp("k54 exists")
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
    bar(timePlot_54,counts_54)
    
    N_tot = sum(counts_54);
    temp_tot = 0;
    surv_prob_54 = [];
    
    for i = 1:numel(counts_54)
        
        temp = counts_54(i);
        temp_tot = temp_tot + temp;
        
        surv_prob_54 = vertcat(surv_prob_54, (1 - (temp_tot / N_tot)));
    end
    
    surv_prob_54_norm = surv_prob_54 / max(surv_prob_54);
  %  time = (0:(length(surv_prob_12_norm)-1) )* res;
    
    subplot(1,2,2)
    plot(timePlot_54, surv_prob_54_norm);
    title("Survival probability: 5 --> 4")
    xlabel("Times (s)")
    ylabel("Probability")
    xlim([0,max(timeBarEdges)+0.1])
    hold off
    
    timeVector_54 = transpose(timePlot_54);
    
    save('survProb_54.mat','timeVector_54','surv_prob_54_norm')
    
else
    disp("k54 is zero")
end
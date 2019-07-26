close all
clear all
clc
%%
% AUTHORS: Claire Albrecht (with ref to Brett Israels' code)
% CREATED: July 2019
% MODIFIED:
%
% PURPOSE: Make histogram of dwell times
%
% INPUT: dwell files (from pathDwell program)
%
% OUTPUT: Histogram of dwell times (next will be to create survival
% probability function from that.
%
% NOTE 1: the dwell times are in frames - not actual time
%           multiply by resolution to get time.
% -------------------------------------------------------------------------


% Find dwell files
scanID = '*';
fileKeyWord = [scanID 'dwell.dat'];
fileNames = dir(fileKeyWord);
numFiles = length(fileNames);

res = 0.03;

for molID = 1:numFiles
    
    dwellFileName = fileNames(molID).name;
    
    if fileNames(molID).bytes > 0
        
        dwelltimes = dlmread(dwellFileName);
        dwell_time_frame = dwelltimes(:,3);
        
    end
    
end

dwell_time = dwell_time_frame * res;
% time = length(dwell_time) * res;



N = histcounts(dwell_time,length(dwell_time));
surv_prob = zeros(length(N));

% Make histogram of dwell times
figure
subplot(1,3,1)
dwellHist = histogram(dwell_time,length(dwell_time));
ylabel('Counts');
xlabel('Dwell time (s)');
title('Dwell time histogram');

%dwellHist = histogram(N,length(N));
hold on

%%
N_tot = sum(N);

% Build survival probability plot
temp_tot = 0;

surv_prob = [];

for i = 1:length(N)
    
    temp = N(i);
    temp_tot = temp_tot + temp;
    
    surv_prob = vertcat(surv_prob, 1 - (temp_tot / N_tot));
end


%[Prob, edges, bins] = histcounts(surv_prob);



surv_prob_norm = surv_prob/max(surv_prob);
% time_bins = 0:length(surv_prob_norm)-1;
time_bins = 0:length(surv_prob_norm)-1;
time = time_bins * res;

subplot(1,3,2)
plot(time_bins,surv_prob_norm);
title('Survival probability vs time bins');
ylabel('Survival Probability');
xlabel('Time bins (counts)');
hold on

subplot(1,3,3)
plot(time, surv_prob_norm);
title('Survival probability vs time ');
ylabel('Survival Probability');
xlabel('Time (s)');
hold off

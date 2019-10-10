function [rates_obs] = rates_obs_finder(mtx_trns)

%--------------------------------------------------------------------------
% Identify the rates that exist based on the TDP
%--------------------------------------------------------------------------
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

%%
%--------------------------------------------------------------------------
% Which rates exist in our data set?
%--------------------------------------------------------------------------
% Define a column of the possible rates
% rates_possible = [k12; k13; k14; k15; k21; k23; k24; k25; k31; k32; k34; k35; k41; k42; k43; k45; k51; k52; k53; k54];

% Ask the data set if the FRET transitions corresponding to each rate exist
%  in the set of rounded data.
%%
rates_obs = {};


% Does rate k12 exist? If so, write k12, if not k12 = 0.
rates_obs{1,1} = 'k12';
if ismember(k12_FRET_a,round_mtx_trns,'rows')||ismember(k12_FRET_b,round_mtx_trns,'rows')
    % k12;
    rates_obs{1,2} = k12;
    disp('k12 exists')
else
    %k12 = 0;
    rates_obs{1,2} = double(0);
    disp('k12 = 0')
end




% Does rate k13 exist? If so, write k13, if not k13 = 0.
rates_obs{2,1} = 'k13';
if ismember(k13_FRET_a,round_mtx_trns,'rows')||ismember(k13_FRET_b,round_mtx_trns,'rows')
    % k13;
    rates_obs{2,2} = k13;
    disp('k13 exists')
else
    % k13 = 0;
    rates_obs{2,2} = double(0);
    disp('k13 = 0')
end



% Does rate k14 exist? If so, write k14, if not k14 = 0.
rates_obs{3,1} = 'k14';
if ismember(k14_FRET_a,round_mtx_trns,'rows')||ismember(k14_FRET_b,round_mtx_trns,'rows')
    %k14;
    rates_obs{3,2} = k14;
    disp('k14 exists')
else
    %k14 = 0;
    rates_obs{3,2} = double(0);
    disp('k14 = 0')
end


% Does rate k15 exist? If so, write k15, if not k15 = 0.
rates_obs{4,1} = 'k15';
if ismember(k15_FRET,round_mtx_trns,'rows')
    % k15;
    rates_obs{4,2} = k15;
    disp('k15 exists')
else
    % k15 = 0;
    rates_obs{4,2} = double(0);
    disp('k15 = 0')
end


% Does rate k21 exist? If so, write k21, if not k21 = 0.
rates_obs{5,1} = 'k21';
if ismember(k21_FRET_a,round_mtx_trns,'rows')||ismember(k21_FRET_b,round_mtx_trns,'rows')
    % k21;
    rates_obs{5,2} = k21;
    disp('k21 exists')
else
    %  k21 = 0;
    rates_obs{5,2} = double(0);
    disp('k21 = 0')
end


% Does rate k23 exist? If so, write k23, if not k23 = 0.
rates_obs{6,1} = 'k23';
if ismember(k23_FRET_aa,round_mtx_trns,'rows')||ismember(k23_FRET_ab,round_mtx_trns,'rows')||ismember(k23_FRET_ba,round_mtx_trns,'rows')||ismember(k23_FRET_bb,round_mtx_trns,'rows')
    % k23;
    rates_obs{6,2} = k23;
    disp('k23 exists')
else
    % k23 = 0;
    rates_obs{6,2} = double(0);
    disp('k23 = 0')
end


% Does rate k24 exist? If so, write k24, if not k24 = 0.
rates_obs{7,1} = 'k24';
if ismember(k24_FRET_aa,round_mtx_trns,'rows')||ismember(k24_FRET_ab,round_mtx_trns,'rows')||ismember(k24_FRET_ba,round_mtx_trns,'rows')||ismember(k24_FRET_bb,round_mtx_trns,'rows')
    % k24;
    rates_obs{7,2} = k24;
    disp('k24 exists')
else
    % k24 = 0;
    rates_obs{7,2} = double(0);
    disp('k24 = 0')
end



% Does rate k25 exist? If so, write k25, if not k25 = 0.
rates_obs{8,1} = 'k25';
if ismember(k25_FRET_a,round_mtx_trns,'rows')||ismember(k25_FRET_b,round_mtx_trns,'rows')
    % k25;
    rates_obs{8,2} = k25;
    disp('k25 exists')
else
    %k25 = 0;
    rates_obs{8,2} = double(0);
    disp('k25 = 0')
end


% Does rate k31 exist? If so, write k31, if not k31 = 0.
rates_obs{9,1} = 'k31';
if ismember(k31_FRET_a,round_mtx_trns,'rows')||ismember(k31_FRET_b,round_mtx_trns,'rows')
    % k31;
    rates_obs{9,2} = k31;
    disp('k31 exists')
else
    %  k31 = 0;
    rates_obs{9,2} = double(0);
    disp('k31 = 0')
end

% if exist('rates')
%     rates = vertcat(rates, [k31]);
% else
%     rates = [k31];
% end



% Does rate k32 exist? If so, write k32, if not k32 = 0.
rates_obs{10,1} = 'k32';
if ismember(k32_FRET_aa,round_mtx_trns,'rows')||ismember(k32_FRET_ab,round_mtx_trns,'rows')||ismember(k32_FRET_ba,round_mtx_trns,'rows')||ismember(k32_FRET_bb,round_mtx_trns,'rows')
  %  k32;
  rates_obs{10,2} = k32;
    disp('k32 exists')
else
   % k32 = 0;
    rates_obs{10,2} = double(0);
    disp('k32 = 0')
end



% Does rate k34 exist? If so, write k34, if not k34 = 0.
rates_obs{11,1} = 'k34';
if ismember(k34_FRET_aa,round_mtx_trns,'rows')||ismember(k34_FRET_ab,round_mtx_trns,'rows')||ismember(k34_FRET_ba,round_mtx_trns,'rows')||ismember(k34_FRET_bb,round_mtx_trns,'rows')
  rates_obs{11,2} = k34;
    disp('k34 exists')
else
    rates_obs{11,2} = double(0);
    disp('k34 = 0')
end


% Does rate k35 exist? If so, write k35, if not k35 = 0.
rates_obs{12,1} = 'k35';
if ismember(k35_FRET_a,round_mtx_trns,'rows')||ismember(k35_FRET_b,round_mtx_trns,'rows')
   rates_obs{12,2} = k35;
    disp('k35 exists')
else
   rates_obs{12,2} = double(0);
    disp('k35 = 0')
end


% Does rate k41 exist? If so, write k41, if not k41 = 0.
rates_obs{13,1} = 'k41';
if ismember(k41_FRET_a,round_mtx_trns,'rows')||ismember(k41_FRET_b,round_mtx_trns,'rows')
   % k41;
   rates_obs{13,2} = k41;
    disp('k41 exists')
else
  %  k41 = 0;
  rates_obs{13,2} = double(0);
    disp('k41 = 0')
end


% Does rate k42 exist? If so, write k42, if not k42 = 0.
rates_obs{14,1} = 'k42';
if ismember(k42_FRET_aa,round_mtx_trns,'rows')||ismember(k42_FRET_ab,round_mtx_trns,'rows')||ismember(k42_FRET_ba,round_mtx_trns,'rows')||ismember(k42_FRET_bb,round_mtx_trns,'rows')
   % k42;
   rates_obs{14,2} = k42;
    disp('k42 exists')
else
   % k42 = 0;
    rates_obs{14,2} = double(0);
    disp('k42 = 0')
end


% Does rate k43 exist? If so, write k43, if not k43 = 0.
rates_obs{15,1} = 'k43';
if ismember(k43_FRET_aa,round_mtx_trns,'rows')||ismember(k43_FRET_ab,round_mtx_trns,'rows')||ismember(k43_FRET_ba,round_mtx_trns,'rows')||ismember(k43_FRET_bb,round_mtx_trns,'rows')
   % k43;
   rates_obs{15,2} = k43;
    disp('k43 exists')
else
   % k43 = 0;
    rates_obs{15,2} = double(0);
    disp('k43 = 0')
end


% Does rate k45 exist? If so, write k45, if not k45 = 0.
rates_obs{16,1} = 'k45';
if ismember(k45_FRET_a,round_mtx_trns,'rows')||ismember(k45_FRET_b,round_mtx_trns,'rows')
   % k45;
    rates_obs{16,2} = k45;
    disp('k45 exists')
else
   % k45 = 0;
   rates_obs{16,2} = double(0);
    disp('k45 = 0')
end


% Does rate k51 exist? If so, write k51, if not k51 = 0.
rates_obs{17,1} = 'k51';
if ismember(k51_FRET,round_mtx_trns,'rows')
   % k51;
   rates_obs{17,2} = k51;
    disp('k51 exists')
else
    % k51 = 0;
    rates_obs{17,2} = double(0);
    disp('k51 = 0')
end


% Does rate k52 exist? If so, write k52, if not k52 = 0.
rates_obs{18,1} = 'k52';
if ismember(k52_FRET_a,round_mtx_trns,'rows')||ismember(k52_FRET_b,round_mtx_trns,'rows')
   % k52;
   rates_obs{18,2} = k52;
    disp('k52 exists')
else
   % k52 = 0;
   rates_obs{18,2} = double(0);
    disp('k52 = 0')
end


% Does rate k53 exist? If so, write k53, if not k53 = 0.
rates_obs{19,1} = 'k53';
if ismember(k53_FRET_a,round_mtx_trns,'rows')||ismember(k53_FRET_b,round_mtx_trns,'rows')
   % k53;
   rates_obs{19,2} = k53;
    disp('k53 exists')
else
   % k53 = 0;
   rates_obs{19,2} = double(0);
    disp('k53 = 0')
end


% Does rate k54 exist? If so, write k54, if not k54 = 0.
rates_obs{20,1} = 'k54';
if ismember(k54_FRET_a,round_mtx_trns,'rows')||ismember(k54_FRET_b,round_mtx_trns,'rows')
   % k54;
   rates_obs{20,2} = k54;
    disp('k54 exists')
else
    % k54 = 0;
    rates_obs{20,2} = double(0);
    disp('k54 = 0')
end


%%
%--------------------------------------------------------------------------
% Collect observed rates
%--------------------------------------------------------------------------



[k12, k13, k14, k15, k21, k23, k24, k25, k31, k32, k34, k35, k41, k42, k43, k45, k51, k52, k53, k54] = rates_obs{1:20,2}



%%
save('dwellTDP_trim','rates_obs')


% Note: this ignores the 0.33 and 0.42 transition which I cannot currently
% make sense of. They do not correspond to FRET states and they are very
% small peaks. They may be a byproduct of false transitions from the
% stitching together of files. Still need to fix that.

function guess = GenAlg_v14_1p0uMgp32
% Fitting to mean subtracted 2-pt. and 4-pt. TCFs.

% Fitting real expt.
initial_pop = 400;
generations = 1000;
maxrepeats = 7;

% Number of individuals who reproduce
numreproduce = round(initial_pop/5);

% Number of mutations per generation
mutations = round(initial_pop/5);

% Bounds
t01_bounds = [1,200];
t10_bounds = [1,2000];
t20_bounds = [1,200000];
t12_bounds = [1,2000];
t21_bounds = [1,2000];
t23_bounds = [1,200];
t32_bounds = [1,2000];
val0_bounds = [.80999,.81001];
val1_bounds = [.56,.81];
%val2_bounds = [.56,.56001];
val3_bounds = [.55999,.56001];

% Used 10 to 5000 all times, 0 to 1 all states.

% Maximum mutation value
t01_mutate = 50;
t10_mutate = 400;
t12_mutate = 40000;
t21_mutate = 400;
t20_mutate = 400;
t23_mutate = 50;
t32_mutate = 400;
val0_mutate = 0;
val1_mutate = .025;
%val2_mutate = 0;
val3_mutate = 0;

repeated = 0;

%%%%%%%%%%%%%%%%%%% Parameters for computing rms %%%%%%%%%%%%%%%%%%%%%%%%%%
rmscalc = @multigoaltcf_analytical;

% Weights
w_2pt = 10*3*10000;
% w_4pt_t0 = 5000;
w_4pt_t1 = 10*500;
w_4pt_t5 = 5*500;
w_4pt_t10 = 5*500;
w_4pt_t15 = 5*500;
w_4pt_t20 = 5*500;
w_4pt_t30 = 5*500;
% w_4pt_t40 = 0;%1000;
% w_4pt_t50 = 0;%1000;
% w_4pt_t75 = 0;
% w_4pt_t100 = 0;
w_frethist = 10*3*300; % Weighting for FRET hist comparison
w_eig1 = 30*100;
w_eig2 = 30*100;

% tau1range = [1:100,200:100:500]';
% tau1range = [1:500]';
tau1range = [1:250]';
tau3range = tau1range;

tcftime = [1:.1:500]';

% Weighting functions
% 2-pt. TCF weighting function
% w1func = ones(length(tcftime),1);
w1func = 1./(sqrt(tcftime));
% 4-pt. TCF weighting function
w2func = 1./sqrt(tau1range)*(1./(sqrt(tau1range')));

% Open 0.1 uM gp32 experimental files
% expt2pt = textread('TCFavgNorm_100usRes_19Mol.dat');
% expt2pt = textread('TCFavgNorm_100usRes_1p0uM_DoubleFit.dat');
% expt2pt = expt2pt(:,2);
% expt2pt = expt2pt(1:4991,2);
expt2pt = textread('TCFavgNorm_meansub_DoubleFitNew_1p0uMgp32.dat');
expt2pt = expt2pt/expt2pt(1);
% A = textread('FourCorrFinal_Fit_1p0uMgp32_tau0.dat');
% A = A(tau1range,tau3range);
B = textread('FourCorrFinal_Fit_1p0uMgp32_tau1.dat');
B = B(tau1range,tau3range);
C = textread('FourCorrFinal_Fit_1p0uMgp32_tau5.dat');
C = C(tau1range,tau3range);
C1 = textread('FourCorrFinal_Fit_1p0uMgp32_tau10.dat');
C1 = C1(tau1range,tau3range);
C2 = textread('FourCorrFinal_Fit_1p0uMgp32_tau15.dat');
C2 = C2(tau1range,tau3range);
D = textread('FourCorrFinal_Fit_1p0uMgp32_tau20.dat');
D = D(tau1range,tau3range);
D1 = textread('FourCorrFinal_Fit_1p0uMgp32_tau30.dat');
D1 = D1(tau1range,tau3range);
% D2 = textread('FourCorrFinal_1p0uMgp32_19mol_tau40.dat');
% D2 = D2(tau1range,tau3range);
% E = textread('FourCorrFinal_1p0uMgp32_19mol_tau50.dat');
% E = E(tau1range,tau3range);
% E1 = textread('FourCorrFinal_1p0uMgp32_19mol_tau75.dat');
% E1 = E1(tau1range,tau3range);
% F = textread('FourCorrFinal_1p0uMgp32_19mol_tau100.dat');
% F = F(tau1range,tau3range);

% A = A/A(1,1);
B = B/B(1,1);
C = C/C(1,1);
C1 = C1/C1(1,1);
C2 = C2/C2(1,1);
D = D/D(1,1);
D1 = D1/D1(1,1);
% D2 = D2/D2(1,1);
% E = E/E(1,1);
% E1 = E1/E1(1,1);
% F = F/F(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End rms parameters %%%%%%%%%%%%%%%%%%%%

% Generation 1
population = rand(initial_pop,10);
population(:,1) = population(:,1)*(t01_bounds(2) - t01_bounds(1)) + t01_bounds(1);
population(:,2) = population(:,2)*(t10_bounds(2) - t10_bounds(1)) + t10_bounds(1);
population(:,3) = population(:,3)*(t20_bounds(2) - t20_bounds(1)) + t20_bounds(1);
population(:,4) = population(:,4)*(t12_bounds(2) - t12_bounds(1)) + t12_bounds(1);
population(:,5) = population(:,5)*(t21_bounds(2) - t21_bounds(1)) + t21_bounds(1);
population(:,6) = population(:,6)*(t23_bounds(2) - t23_bounds(1)) + t23_bounds(1);
population(:,7) = population(:,7)*(t32_bounds(2) - t32_bounds(1)) + t32_bounds(1);
population(:,8) = population(:,8)*(val0_bounds(2) - val0_bounds(1)) + val0_bounds(1);
population(:,9) = population(:,9)*(val1_bounds(2) - val1_bounds(1)) + val1_bounds(1);
population(:,10) = population(:,10)*(val3_bounds(2) - val3_bounds(1)) + val3_bounds(1);

newpopulation = zeros(size(population));

popfit = zeros(initial_pop,1);

guess = zeros(generations,11);

tic
parfor k = 1:initial_pop
   popfit(k) = rmscalc(population(k,8),population(k,9),population(k,10),1/population(k,1),1/population(k,2),1/population(k,3),1/population(k,4),1/population(k,5),1/population(k,6),1/population(k,7));
end
toc

% Evaluate Generation 1
% Sort from highest rated to lowest rated individuals
[popfit,index] = sort(popfit);
population = population(index,:);
guess(1,1:10) = population(1,:); % Current best guess
guess(1,11) = popfit(1); % Current best guess' fit value
population(1,:); % Display best guess to screen
popfit(1)

for m = 2:generations
tic
    % Mix genes of top rated individuals
    for n = 1:initial_pop - 2
        newpopulation(n,:) = [population(ceil(numreproduce*rand(1,1)),1),population(ceil(numreproduce*rand(1,1)),2),population(ceil(numreproduce*rand(1,1)),3),population(ceil(numreproduce*rand(1,1)),4),population(ceil(numreproduce*rand(1,1)),5),population(ceil(numreproduce*rand(1,1)),6),population(ceil(numreproduce*rand(1,1)),7),population(ceil(numreproduce*rand(1,1)),8),population(ceil(numreproduce*rand(1,1)),9),population(ceil(numreproduce*rand(1,1)),10)];
    end
    
    % Keep 2 best parents
    newpopulation(initial_pop - 1,:) = population(2,:);
    newpopulation(initial_pop,:) = population(1,:);
    
    population = newpopulation;
    
    % Mutate some individuals
    for q = 1:mutations
        rnd1 = rand(1,1);
        rnd2 = rand(1,1);
        rnd3 = rand(1,1);
        rnd4 = rand(1,1);
        rnd5 = rand(1,1);
        rnd6 = rand(1,1);
        rnd7 = rand(1,1);
        rnd8 = rand(1,1);
        rnd9 = rand(1,1);
        rnd10 = rand(1,1);
           population(ceil((initial_pop - 2)*rnd1),1) = population(ceil((initial_pop - 2)*rnd1),1) + (2*(rand(1,1) - .5))*t01_mutate;
        if population(ceil((initial_pop - 2)*rnd1),1) > t01_bounds(2)
            population(ceil((initial_pop - 2)*rnd1),1) = t01_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd1),1) < t01_bounds(1)
            population(ceil((initial_pop - 2)*rnd1),1) = t01_bounds(1);
        end
        
        population(ceil((initial_pop - 2)*rnd2),2) = population(ceil((initial_pop - 2)*rnd2),2) + (2*(rand(1,1) - .5))*t10_mutate;
        if population(ceil((initial_pop - 2)*rnd2),2) > t10_bounds(2)
            population(ceil((initial_pop - 2)*rnd2),2) = t10_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd2),2) < t10_bounds(1)
            population(ceil((initial_pop - 2)*rnd2),2) = t10_bounds(1);
        end

        population(ceil((initial_pop - 2)*rnd3),3) = population(ceil((initial_pop - 2)*rnd3),3) + (2*(rand(1,1) - .5))*t20_mutate;
        if population(ceil((initial_pop - 2)*rnd3),3) > t20_bounds(2)
            population(ceil((initial_pop - 2)*rnd3),3) = t20_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd3),3) < t20_bounds(1)
            population(ceil((initial_pop - 2)*rnd3),3) = t20_bounds(1);
        end
        
        population(ceil((initial_pop - 2)*rnd4),4) = population(ceil((initial_pop - 2)*rnd4),4) + (2*(rand(1,1) - .5))*t12_mutate;
        if population(ceil((initial_pop - 2)*rnd4),4) > t12_bounds(2)
            population(ceil((initial_pop - 2)*rnd4),4) = t12_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd4),4) < t12_bounds(1)
            population(ceil((initial_pop - 2)*rnd4),4) = t12_bounds(1);
        end
        
        population(ceil((initial_pop - 2)*rnd5),5) = population(ceil((initial_pop - 2)*rnd5),5) + (2*(rand(1,1) - .5))*t21_mutate;
        if population(ceil((initial_pop - 2)*rnd5),5) > t21_bounds(2)
            population(ceil((initial_pop - 2)*rnd5),5) = t21_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd5),5) < t21_bounds(1)
            population(ceil((initial_pop - 2)*rnd5),5) = t21_bounds(1);
        end
        
        population(ceil((initial_pop - 2)*rnd6),6) = population(ceil((initial_pop - 2)*rnd6),6) + (2*(rand(1,1) - .5))*t23_mutate;
        if population(ceil((initial_pop - 2)*rnd6),6) > t23_bounds(2)
            population(ceil((initial_pop - 2)*rnd6),6) = t23_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd6),6) < t23_bounds(1)
            population(ceil((initial_pop - 2)*rnd6),6) = t23_bounds(1);
        end
        
        population(ceil((initial_pop - 2)*rnd7),7) = population(ceil((initial_pop - 2)*rnd7),7) + (2*(rand(1,1) - .5))*t32_mutate;
        if population(ceil((initial_pop - 2)*rnd7),7) > t32_bounds(2)
            population(ceil((initial_pop - 2)*rnd7),7) = t32_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd7),7) < t32_bounds(1)
            population(ceil((initial_pop - 2)*rnd7),7) = t32_bounds(1);
        end
       
        population(ceil((initial_pop - 2)*rnd8),8) = population(ceil((initial_pop - 2)*rnd8),8) + (2*(rand(1,1) - .5))*val0_mutate;
        if population(ceil((initial_pop - 2)*rnd8),8) > val0_bounds(2)
            population(ceil((initial_pop - 2)*rnd8),8) = val0_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd8),8) < val0_bounds(1)
            population(ceil((initial_pop - 2)*rnd8),8) = val0_bounds(1);
        end
        
        population(ceil((initial_pop - 2)*rnd9),9) = population(ceil((initial_pop - 2)*rnd9),9) + (2*(rand(1,1) - .5))*val1_mutate;
        if population(ceil((initial_pop - 2)*rnd9),9) > val1_bounds(2)
            population(ceil((initial_pop - 2)*rnd9),9) = val1_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd9),9) < val1_bounds(1)
            population(ceil((initial_pop - 2)*rnd9),9) = val1_bounds(1);
        end
        
        population(ceil((initial_pop - 2)*rnd10),10) = population(ceil((initial_pop - 2)*rnd10),10) + (2*(rand(1,1) - .5))*val3_mutate;
        if population(ceil((initial_pop - 2)*rnd10),10) > val3_bounds(2)
            population(ceil((initial_pop - 2)*rnd10),10) = val3_bounds(2);
        elseif population(ceil((initial_pop - 2)*rnd10),10) < val3_bounds(1)
            population(ceil((initial_pop - 2)*rnd10),10) = val3_bounds(1);
        end
        
    end
    
    % Calculate fitness
    parfor k = 1:initial_pop
       popfit(k) = rmscalc(population(k,8),population(k,9),population(k,10),1/population(k,1),1/population(k,2),1/population(k,3),1/population(k,4),1/population(k,5),1/population(k,6),1/population(k,7));
    end

    % Evaluate Generation
    % Sort from highest rated to lowest rated individuals
    [popfit,index] = sort(popfit);
    population = population(index,:);

    guess(m,1:10) = population(1,:); % Current best guess
    guess(m,11) = popfit(1); % Current best guess' fit value
    
    %population(1,:) % Display best guess to screen
    popfit(1)
    
    if (guess(m-1,11) - guess(m,11))/guess(m,11) <= .01
        repeated = repeated + 1;
    else
        repeated = 0;
    end
    if repeated == maxrepeats
        guess = guess(1:m,:);
        break
    end
toc
end



function rms = multigoaltcf_analytical(val0,val1,val3,k01,k10,k20,k12,k21,k23,k32)
    val2 = val1;
rms_array = zeros(15,1);

% Calculate 2-pt. TCF using analytical formula
tcf = TCF_4state(tcftime,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
tcf = tcf/tcf(1);

% Calculate rms deviation from 2-pt. TCF of guess
rms_array(1) = sum(((tcf - expt2pt).*w1func).^2)*w_2pt;

% Calculate 4-pt. TCF from analytical formula with various tau2 times.
% [~,FourCorr0] = FourPtTCF_FourState_v3(tau1range,0,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
% FourCorr0 = FourCorr0/kap0(1,1);
[~,FourCorr1] = FourPtTCF_FourState_v3(tau1range,1,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
%FourCorr1 = FourCorr1/kap1(1,1);
[~,FourCorr5] = FourPtTCF_FourState_v3(tau1range,5,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
%FourCorr5 = FourCorr5/kap5(1,1);
[~,FourCorr10] = FourPtTCF_FourState_v3(tau1range,10,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr15] = FourPtTCF_FourState_v3(tau1range,15,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr20] = FourPtTCF_FourState_v3(tau1range,20,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
%FourCorr20 = FourCorr20/kap20(1,1);
[~,FourCorr30] = FourPtTCF_FourState_v3(tau1range,30,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
% [kap40,FourCorr40] = FourPtTCF_FourState_v3(tau1range,40,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
% [kap50,FourCorr50] = FourPtTCF_FourState_v3(tau1range,50,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
% %FourCorr50 = FourCorr50/kap50(1,1);
% [kap75,FourCorr75] = FourPtTCF_FourState_v3(tau1range,75,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
% [kap100,FourCorr100] = FourPtTCF_FourState_v3(tau1range,100,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
% %FourCorr100 = FourCorr100/kap100(1,1);

% FourCorr0 = FourCorr0/FourCorr0(1,1);
FourCorr1 = FourCorr1/FourCorr1(1,1);
FourCorr5 = FourCorr5/FourCorr5(1,1);
FourCorr10 = FourCorr10/FourCorr10(1,1);
FourCorr15 = FourCorr15/FourCorr15(1,1);
FourCorr20 = FourCorr20/FourCorr20(1,1);
FourCorr30 = FourCorr30/FourCorr30(1,1);
% FourCorr40 = FourCorr40/FourCorr40(1,1);
% FourCorr50 = FourCorr50/FourCorr50(1,1);
% FourCorr75 = FourCorr75/FourCorr75(1,1);
% FourCorr100 = FourCorr100/FourCorr100(1,1);

% Calculate rms deviation from 4-pt. TCF of guess with t2 = 0
% rms_array(2) = sum(sum((A - FourCorr0).^2.*w2func))*w_4pt_t0;
rms_array(3) = sum(sum((B - FourCorr1).^2.*w2func))*w_4pt_t1;
rms_array(4) = sum(sum((C - FourCorr5).^2.*w2func))*w_4pt_t5;
rms_array(5) = sum(sum((D - FourCorr20).^2.*w2func))*w_4pt_t20;
% rms_array(6) = sum(sum((E - FourCorr50).^2.*w2func))*w_4pt_t50;
% rms_array(7) = sum(sum((F - FourCorr100).^2.*w2func))*w_4pt_t100;
rms_array(11) = sum(sum((C1 - FourCorr10).^2.*w2func))*w_4pt_t10;
rms_array(12) = sum(sum((C2 - FourCorr15).^2.*w2func))*w_4pt_t15;
rms_array(13) = sum(sum((D1 - FourCorr30).^2.*w2func))*w_4pt_t30;
% rms_array(14) = sum(sum((D2 - FourCorr40).^2.*w2func))*w_4pt_t40;
% rms_array(15) = sum(sum((E1 - FourCorr75).^2.*w2func))*w_4pt_t75;
rms_array(2) = 0;
rms_array(6) = 0;
rms_array(7) = 0;
rms_array(14) = 0;
rms_array(15) = 0;

% rms_array(2) = sum(sum((A - FourCorr0).^2))*w_4pt_t0;
% rms_array(3) = sum(sum((B - FourCorr1).^2))*w_4pt_t1;
% rms_array(4) = sum(sum((C - FourCorr5).^2))*w_4pt_t5;
% rms_array(5) = sum(sum((D - FourCorr20).^2))*w_4pt_t20;
% rms_array(6) = sum(sum((E - FourCorr50).^2))*w_4pt_t50;
% rms_array(7) = sum(sum((F - FourCorr100).^2))*w_4pt_t100;
% rms_array(11) = sum(sum((C1 - FourCorr10).^2))*w_4pt_t10;
% rms_array(12) = sum(sum((C2 - FourCorr15).^2))*w_4pt_t15;
% rms_array(13) = sum(sum((D1 - FourCorr30).^2))*w_4pt_t30;
% rms_array(14) = sum(sum((D2 - FourCorr40).^2))*w_4pt_t40;
% rms_array(15) = sum(sum((E1 - FourCorr75).^2))*w_4pt_t75;


% Comparing to 1.0 uM FRET histogram
bins = [0:.01:1]';
fithist = 686.9/1271*exp(-((bins-0.8045)/.09628).^2) + exp(-((bins-0.5682)/.1314).^2);
[p0_eq,p1_eq,p2_eq,p3_eq,eig1,eig2,eig3] = FourPtTCF_FourState_xo(tau1range,1,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
val0scaled = 0.79;
val3scaled = 0.587;
xscale = (val1 - val3)/(val0 - val3);
val1scaled = val3scaled + xscale*(val0scaled - val3scaled);
val2scaled = val1scaled;
analhist = p0_eq*exp(-((bins-val0scaled)/.1).^2) + p1_eq*exp(-((bins-val1scaled)/.1).^2) + p2_eq*exp(-((bins-val2scaled)/.1).^2) + p3_eq*exp(-((bins-val3scaled)/.1).^2);
normsy = max(analhist);
analhist = (p0_eq*exp(-((bins-val0scaled)/.1).^2) + p1_eq*exp(-((bins-val1scaled)/.1).^2) + p2_eq*exp(-((bins-val2scaled)/.1).^2) + p3_eq*exp(-((bins-val3scaled)/.1).^2))/normsy;
rms_array(8) = sum(w_frethist*(fithist - analhist).^2);

% Compare to "known" eigenvalues
dif1 = abs(eig1 - 0.072);
dif2 = abs(eig2 - 0.072);
dif3 = abs(eig3 - 0.072);
difa = abs(eig1 - 0.0106);
difb = abs(eig2 - 0.0106);
difc = abs(eig3 - 0.0106);
if dif1 <= dif2
    if dif1 <= dif3
        rms_array(9) = dif1*w_eig1;
        firstis = 1;
        difaa = difb;
        difbb = difc;
    else
        rms_array(9) = dif3*w_eig1;
        firstis = 3;
        difaa = difa;
        difbb = difb;
    end
elseif dif2 <= dif3
    rms_array(9) = dif2*w_eig1;
    firstis = 2;
    difaa = difa;
    difbb = difc;
else
    rms_array(9) = dif3*w_eig1;
    firstis = 3;
    difaa = difa;
    difbb = difb;
end

if difaa <= difbb
    rms_array(10) = difaa*w_eig2;
else
    rms_array(10) = difbb*w_eig2;
end

rms = sum(rms_array);

end

end
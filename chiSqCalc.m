function [chisquared,chisquared_array, chisquared_perGen] = chiSqCalc(rates, A, P, sigma_A, C2_time,C4_time)

global normalizeMode diagnoseMode
global targetHistogram weightingFactor_FREThist FRET_bins
global C2_exp_x C2_exp_y weightingFactor_C2 weightC2func
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
global fitHistMode fitC2Mode fitC4Mode
global sigma_A1 sigma_A2 sigma_A3 sigma_A4 sigma_A5 sigma_A6 sigma_A7 sigma_A8 %sigma_A9
global yoff zoff genNum
global showProgressOnFit_mode

global  K % P
%Initialize an array to hold the chi squared values
chisquared_array = zeros(3,1);

% sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4; sigma_A5; sigma_A6; sigma_A7; sigma_A8];
%------------------------------------------------------------------
% (1) Optimization Target #1: Get a single value for the entire
% hitogram: rms_array(1) (OPTIMIZATION FUNCTION:% multigoaltcf_analytical_3state)
%------------------------------------------------------------------
if fitHistMode == 1
    [Peq,~,hist_sim] = histMaker_Nstate(P, A, sigma_A);
    chisquared_array(1) = sum((reshape(hist_sim,100,1)-targetHistogram).^2)*weightingFactor_FREThist;

    if showProgressOnFit_mode == 1
        % Update figure with guesses to the fit
        figure(1);
        
        hist_sim_temp_plot = plot(FRET_bins,hist_sim,'r-','LineWidth',3);
        if exist('hist_sim_temp_plot','var') == 1
            pause(0.01);
            delete(hist_sim_temp_plot);
        end
        
    end
end

%------------------------------------------------------------------------------------------
% (2) Optimization Target #2: Get a single value for the entire 2-point TCF: rms_array(1)
%-----------------------------------------------------------------------------------------
if fitC2Mode == 1
    time = C2_time;
    [P_c2, ~, ~, time] = k2P(K,time);
    [time, C2_sim, ~] = P2C(P_c2, K, time);
    C2_sim = C2_sim + yoff;
    if normalizeMode == 1
        C2_sim = C2_sim./C2_sim(1);
        end
    chisquared_array(2) = sum(((reshape(C2_sim,length(C2_sim),1) - C2_exp_y).*weightC2func).^2)*weightingFactor_C2;
    time = time_sim;
    if showProgressOnFit_mode == 1
        % Update figure with guesses to the fit
        figure(2);
        
        C2_sim_temp_plot = plot(C2_exp_x,C2_sim,'r-','LineWidth',3);
        if exist('C2_shift_temp_plot','var') == 1
            pause(0.01);
            delete(C2_shift_temp_plot);
        end
        
    end
end

%--------------------------------------------------------------------------
% Optimization Target #3: Get a single value for the entire 4-point TCF: rms_array(3)
%--------------------------------------------------------------------------
if fitC4Mode == 1
    tau2 = 0;
    time = C4_time;
    [P_c4, ~, ~, time] = k2P(K,time);
    [time, ~, C4_sim] = P2C(P_c4, K, time);
    C4_tau2eq0_sim = C4_sim + zoff;
    if normalizeMode == 1
        C4_tau2eq0_sim = C4_tau2eq0_sim./C4_tau2eq0_sim(1,1);
    end
    chisquared_array(3) = mean(mean(((C4_tau2eq0_sim - C4_tau2eq0_exp).^2).*wC4func))*weightingFactor_C4_t0;
    time = time_sim;
    if showProgressOnFit_mode == 1
        % Update figure with guesses to the fit
        figure(3);
        
        C4_tau2eq0_sim_temp_plot = surf(C4_tau1range,C4_tau1range,C4_tau2eq0_sim);
        if exist('C4_tau2eq0_sim_temp_plot','var') == 1
            pause(0.01);
            delete(C4_tau2eq0_sim_temp_plot);
            
        end
        
    end
end
    
%--------------------------------------------------------------------------
% Sum all of the indivual root-mean-squares to get an overall rms
%--------------------------------------------------------------------------
chisquared = sum(chisquared_array);

if genNum  == 1
    chisquared_perGen  = [];
    chisquared_perGen(genNum) = chisquared;
else
    chisquared_perGen(genNum) = chisquared;
end

% if genNum == 1
% chisquared_perGen = [];
% chisquared = sum(chisquared_array);
% chisquared_perGen =  [chisquared_perGen; chisquared];
% 
% else
% chisquared = sum(chisquared_array);
% chisquared_perGen =  [chisquared_perGen; chisquared];
% end

% disp('size of chisquared_array')
% disp(size(chisquared_array))
%------------------------------------------------------------------
% Display relative contributions
%------------------------------------------------------------------
if diagnoseMode == 1
    %Plot the relative contribution of each Chi-squared paramater
    %     histContribution = 0;
    %     C2Contriubtion = 0;
    %     C4Contribution = 0;
    if fitHistMode == 1
        histContribution = chisquared_array(1)*100/chisquared;
        disp(['     The contribution of the Histogram to the rms_array is ' num2str(histContribution) '%']);
    end
    if fitC2Mode == 1
        C2Contriubtion = chisquared_array(2)*100/chisquared;
        disp(['     The contribution of the 2ptTCF to the rms_array is ' num2str(C2Contriubtion) '%']);
    end
    if fitC4Mode == 1
        C4Contribution = chisquared_array(3)*100/chisquared;
        disp(['     The contribution of the 4-point TCF to the rms_array is ' num2str(C4Contribution) '%']);
    end
    
    
end

end
    
    
    

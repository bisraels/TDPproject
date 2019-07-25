%__________________________________________________________________________
% Author: Brett Israels & Claire
%
% FUNCTION: hmmRateFinder
%
% PURPOSE: Takes the equilibrium populations from the HMM and returns the
% rates which could have given rise to it based on solns to master eqn.
%
% INPUT: (1) Solutions to the master equation (Sym_prob.mat);
%
% OUTPUT: (1) Folder named after the function
%         (2) Plot of Probability as a function of rates
%
% DEPENDENCIES: (1) Needs the Solutions to the master equation
%               (2) Needs the HMM fit (path-file)
%
% ALGORITHM:
%  1) Run HMM and Determine the numer of states and the state-to-state
%  connections (A := Adjacency Matrix)
% 2) Use the HMM fit to find P_i^{eq} for each state (Target function)
%    (Use ____.m --> exp_prob.mat)
% 3) Solve the rate matrix master equation for P (the equilibrium state
%  populations)
%       (TDP_probcalc.m)
% 4) Make initial guesses for the rates (unknown quantitites)
% 5) Use the solution and the current guess parameters to calculate Pi_eq
%       (TDP_probcalc_ratevals.m)
% 6) Calculate chi_square (Diffrence between Pi_exp and Pi_calc)
%           in function "minimizeThis"
% 7) Change values of kij to minimize chi_square until chi_square <
% threshold (0.001)
%
% MODIFICATION LOG:
% BI 20190722 Creation of fitter.
% BI 20190723 Initializing global options
%
%
%__________________________________________________________________________
functionName = 'hmmRateFinder';

%--------------------------------------------------------------------------
% Navigate to the folder with the data
%--------------------------------------------------------------------------
cd('/Users/bisraels/Dropbox/programs/SingleMolecule/TDP/'); %Folder containting all relevant files

%--------------------------------------------------------------------------
% Program Paramameters
%--------------------------------------------------------------------------
global Nstates EqPopArray_exp %EqPopArray_exp is a structure. Access attributes with dot operator.
Nstates = 5; %The number of states found by HMM
saveMode = 1;
saveFigMode = 1;
savePngMode = 1;

add_date_str_mode = 0;
special_folder_ending = 'fmincon';
%--------------------------------------------------------------------------
% Prepare for program output
%--------------------------------------------------------------------------
if saveMode
    %Optional: Add a date-string to end of file
    if add_date_str_mode
        date_str = ['_' datestr(now,'yyyymmdd')];
    else
        date_str = '';
    end
    
    %OPTIONAL: Add a special folder ending
    special_folder_ending = ['_' special_folder_ending];
    if length(special_folder_ending) < 2
        special_folder_ending = '';
    end
    
    %CREATE the folder if its not there
    outputFolderName = [functionName '_output' date_str special_folder_ending];
    %If the folder doesn't exist, create it
    if exist(outputFolderName,'dir') ~= 7
        fprintf('     ***Making a folder called %s to store output.\r\n',outputFolderName);
        mkdir(outputFolderName);
    else
        disp('Already detected an outout folder');
    end
end

%--------------------------------------------------------------------------
% Load the target function (Results of HMM Fitting)
%--------------------------------------------------------------------------
hmmResultsFilename = 'exp_prob.mat';
EqPopArray_exp = load(hmmResultsFilename,'FRET_states','EqPopArray');  %EqPopArray is the equilibrium populations sorted by FRET state (low --> High)

%--------------------------------------------------------------------------
% Plot the results from the HMM
%--------------------------------------------------------------------------
figure(1)
clf;
fig_1 = bar(EqPopArray_exp.FRET_states,EqPopArray_exp.EqPopArray);
fig_1.FaceAlpha = .5;
xlabel('FRET Value');
ylabel('Probability');
set(gcf,'Color','w');
legend({'HMM FRET'})
hold on;
drawnow();
%--------------------------------------------------------------------------
%  Set initial guesses for the parameters (x0)
%--------------------------------------------------------------------------
if Nstates == 5
    x0 = 1./randi([1,1000],1,Nstates*(Nstates-1))*1e-3;
end

%**************************************************************************
% Perform the minimization of chisquare using the array of parameters
%**************************************************************************
Ntrials = 1;
figure(2)
for i = 1:Ntrials

options = optimset('Display','iter','PlotFcns',@optimplotfval);
% [x,fval] = fminsearch(@minimizeThis,x0,options);
% % x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
[x,fval] = fmincon(@minimizeThis,x0,[],[],[],[],zeros(size(x0)),[],[],options);
fvalNew = fval;
%**************************************************************************
% Save the results
%**************************************************************************
if saveMode
    
    foutName = ['hmmRateFinder_Results'];
    filepath_fit = [outputFolderName filesep() foutName '.mat'];
    
    if exist(filepath_fit,'file') ~= 2
        disp('Did not find the file. Saving the initial fit');
        iterNumber = 1;
        save(filepath_fit,'fval','iterNumber','x');
    else%The file already exist
        load(filepath_fit,'fval','iterNumber');
        iterNumber = iterNumber + 1 ;
        %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
        if fvalNew < fval
            disp(['*** Found a better fit: ' num2str(fvalNew) ' < ' num2str(fval) '. Saving.']);
            fval = fvalNew;
            save(filepath_fit,'fval','iterNumber','x','iterNumber');
        else
            %                         disp(['     appending the iteration number to' filepath]);
            save(filepath_fit,'iterNumber','-append')
        end
    end
    
end

%**************************************************************************
% Print the Results
%**************************************************************************
switch  Nstates
    case 5
        k12 = x(1);
        k13 = x(2);
        k14 = x(3);
        k15 = x(4);
        k21 = x(5);
        k23 = x(6);
        k24 = x(7);
        k25 = x(8);
        k31 = x(9);
        k32 = x(10);
        k34 = x(11);
        k35 = x(12);
        k41 = x(13);
        k42 = x(14);
        k43 = x(15);
        k45 = x(16);
        k51 = x(17);
        k52 = x(18);
        k53 = x(19);
        k54 = x(20);
        masterEqnSolutions_filename = 'sym_prob.mat';
        load(masterEqnSolutions_filename,'P1_n5','P2_n5','P3_n5','P4_n5','P5_n5')
        eqPopArray_calc(1) = vpa(subs(P1_n5));
        eqPopArray_calc(2) = vpa(subs(P2_n5));
        eqPopArray_calc(3) = vpa(subs(P3_n5));
        eqPopArray_calc(4) = vpa(subs(P4_n5));
        eqPopArray_calc(5) = vpa(subs(P5_n5));
        
        
        
end

%Plot the calculated equilibrium population array
figure(1)
hold on;
fig_2 = bar(EqPopArray_exp.FRET_states,eqPopArray_calc);
fig_2.FaceColor = 'r';
fig_2.FaceAlpha = .5;
fig_2.BarWidth = .5;

xlabel('FRET Value');
ylabel('Probability');
set(gcf,'Color','w');
legend({'HMM FRET','P(k_{ij})'})
hold off;
drawnow();
end
if saveFigMode
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    foutName =  ['ProbHistorgram'];
    saveas(gcf,[outputFolderName filesep() foutName],['fig']);
    disp(['     Saving figure as ' foutName '.fig in ' outputFolderName]);
    if savePngMode
        extension = 'png';
        print(fig,[outputFolderName filesep() foutName],['-d' extension],'-r0');
        disp(['     Saving figure as ' foutName '.' extension ' in ' outputFolderName]);
    end
end
%**************************************************************************
% Minimzation Function
%**************************************************************************
function [chisquare] = minimizeThis(x0)%x0 is the initial guess for rates

global EqPopArray_exp Nstates

eqPopArray_calc = rate_to_pop(x0);

%Function to be minimized
chisquare = sum((EqPopArray_exp.EqPopArray - eqPopArray_calc).^2);
end


%**************************************************************************
% Calculate the equilibrium population from solutions to master equation
%**************************************************************************
function eqPopArray_calc = rate_to_pop(x0)
masterEqnSolutions_filename = 'sym_prob.mat';
global Nstates
eqPopArray_calc = zeros(1,Nstates);

switch Nstates
    case 2
        k12 = x0(1);
        k21 = x0(2);
        
        load(masterEqnSolutions_filename,'P1_n2','P2_n2')
        eqPopArray_calc(1) = vpa(subs(P1_n2));
        eqPopArray_calc(2) = vpa(subs(P2_n2));
        
    case 3
        k12 = x0(1);
        k13 = x0(2);
        k21 = x0(3);
        k23 = x0(4);
        k31 = x0(5);
        k32 = x0(6);
        
        load(masterEqnSolutions_filename,'P1_n3','P2_n3','P3_n3')
        eqPopArray_calc(1) = vpa(subs(P1_n3));
        eqPopArray_calc(2) = vpa(subs(P2_n3));
        eqPopArray_calc(3) = vpa(subs(P3_n3));
        
    case  5
        k12 = x0(1);
        k13 = x0(2);
        k14 = x0(3);
        k15 = x0(4);
        k21 = x0(5);
        k23 = x0(6);
        k24 = x0(7);
        k25 = x0(8);
        k31 = x0(9);
        k32 = x0(10);
        k34 = x0(11);
        k35 = x0(12);
        k41 = x0(13);
        k42 = x0(14);
        k43 = x0(15);
        k45 = x0(16);
        k51 = x0(17);
        k52 = x0(18);
        k53 = x0(19);
        k54 = x0(20);
        
        load(masterEqnSolutions_filename,'P1_n5','P2_n5','P3_n5','P4_n5','P5_n5')
        eqPopArray_calc(1) = vpa(subs(P1_n5));
        eqPopArray_calc(2) = vpa(subs(P2_n5));
        eqPopArray_calc(3) = vpa(subs(P3_n5));
        eqPopArray_calc(4) = vpa(subs(P4_n5));
        eqPopArray_calc(5) = vpa(subs(P5_n5));
end


end
% function [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix] = fileFinder(computer_terminal_str, terminalID, constructName, modelName,protein_str)
function [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix] = fileFinder(terminalID, constructName, modelName, protein_str)

global verboseMode

% PURPOSE:  (1) finds GenAlg files 
%           (2) loads targetdata
%
% INPUT:    terminalID from computerMode (don't need computer_terminal_str)
%           constructName, modelName, protein_str (in GenAlg script)

% verboseMode = 1;

% Part 1: Find the folder for the GenAlg you are running
%--------------------------------------------------------------------------
% Navigate to the correct folder
%--------------------------------------------------------------------------
if contains(terminalID, 'C:\')          % Windows computer
    if verboseMode == 1
        disp('      Entering windows mode...')
    end
    % PART 1: Find the folder for the GenAlg you are running
    GenAlgFilesLocationPrefix =  [terminalID '\Dropbox\MarcusLab\Data\smData_Processed\' constructName '\' protein_str '\ChosenMolecules\genAlgFiles\'  modelName];
    
    % PART 2: Load the target data
    genAlgDataLocationPrefix =  [terminalID 'Dropbox\MarcusLab\Data\smData_Processed\' constructName '\' protein_str '\ChosenMolecules\genAlgFiles\TargetData'];

elseif contains(terminalID,'/Users/')   % Mac computer
%     if verboseMode == 1
%         disp('      Entering Mac mode...')
%     end
    % PART 1: Find the folder for the GenAlg you are running
    GenAlgFilesLocationPrefix =  [terminalID '/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/' protein_str '/ChosenMolecules/genAlgFiles/' modelName];
    
    % PART 2: Load the target data
    genAlgDataLocationPrefix =  [terminalID 'Dropbox/MarcusLab/Data/smData_Processed/' constructName '/' protein_str '/ChosenMolecules/genAlgFiles/TargetData'];

else
    error('Cannot detect computer')
end

if exist(GenAlgFilesLocationPrefix,'dir') ~= 7
    disp('Making a directory to hold the outputs of the model');
    mkdir(GenAlgFilesLocationPrefix);
end
cd(GenAlgFilesLocationPrefix);
wd = pwd;
if verboseMode == 1
    disp(['You are now in the folder: ' wd]);
end


if exist(genAlgDataLocationPrefix,'dir') ~= 7
    error(['Cannot locate data folder ' genAlgDataLocationPrefix]);
end

if verboseMode == 1
    disp(['Looking for data in ' genAlgDataLocationPrefix]);
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % If you want the two parts separate:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Part 1: Find the folder for the GenAlg you are running
% %--------------------------------------------------------------------------
% % Navigate to the correct folder
% %--------------------------------------------------------------------------
% if contains(terminalID, 'C:\')          % Windows computer
%     if verboseMode == 1
%         disp('      Entering windows mode...')
%     end
%     % PART 1: Find the folder for the GenAlg you are running
%     GenAlgFilesLocationPrefix =  [terminalID '\Dropbox\MarcusLab\Data\smData_Processed\' constructName '\' protein_str '\ChosenMolecules\genAlgFiles\'  modelName];
% 
% elseif contains(terminalID,'/Users/')   % Mac computer
%     if verboseMode == 1
%         disp('      Entering Mac mode...')
%     end
%     % PART 1: Find the folder for the GenAlg you are running
%     GenAlgFilesLocationPrefix =  [terminalID '/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/' protein_str '/ChosenMolecules/genAlgFiles/' modelName];
% end
% 
% if exist(GenAlgFilesLocationPrefix,'dir') ~= 7
%     disp('Making a directory to hold the outputs of the model');
%     mkdir(GenAlgFilesLocationPrefix);
% end
% cd(GenAlgFilesLocationPrefix);
% wd = pwd;
% if verboseMode == 1
%     disp(['You are now in the folder: ' wd]);
% end
% 
% 
% % PART 2: Load the target data (and plot if opted)
% %--------------------------------------------------------------------------
% % Filenames of the data (PART 2: Load the target data)
% %--------------------------------------------------------------------------
% 
% if contains(terminalID, 'C:\')          % Windows
%     genAlgDataLocationPrefix =  [terminalID 'Dropbox\MarcusLab\Data\smData_Processed\' constructName '\' protein_str '\ChosenMolecules\genAlgFiles\targetData'];
% elseif contains(terminalID,'/Users/')   % Mac
%     genAlgDataLocationPrefix =  [terminalID 'Dropbox/MarcusLab/Data/smData_Processed/' constructName '/' protein_str '/ChosenMolecules/genAlgFiles/targetData'];
% else
%     error('Cannot detect computer')
% end
% 
% if exist(genAlgDataLocationPrefix,'dir') ~= 7
%     error(['Cannot locate data folder ' genAlgDataLocationPrefix]);
% end
% 
% if verboseMode == 1
%     disp(['Looking for data in ' genAlgDataLocationPrefix]);
% end
% 
% 
% 
% 
% 


%--------------------------------------------------------------------------
% Navigate to the correct folder
%--------------------------------------------------------------------------
% switch computer_terminal_str
%     case 'computer_baimi_mode' %Work Desktop
%         GenAlgFilesLocationPrefix =  ['C:\Users\baimi\Dropbox\MarcusLab\Data\smData_Processed\' constructName '\gp32_0p0uM\ChosenMolecules\genAlgFiles\3state123_cyclical'];
%     case 'computer_bisraels_mode' %Macbook pro
%         GenAlgFilesLocationPrefix =  ['/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
%     case 'computer_claire_mode' %Macbook pro
%         GenAlgFilesLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
% end
% if exist(GenAlgFilesLocationPrefix,'dir') ~= 7
%     disp('Making a directory to hold the outputs of the model');
%     mkdir(GenAlgFilesLocationPrefix);
% end
% cd(GenAlgFilesLocationPrefix);
% wd = pwd;
% if verboseMode == 1
%     disp(['You are now in the folder: ' wd]);
% end




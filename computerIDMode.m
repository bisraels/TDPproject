function [terminalID] = computerIDMode()
% computer_terminal_str is the string used to identify the computer
switch nargin
    case 0
        wd = pwd;
end


% verboseMode = 1;

% Define identifying strings in each computer's command line
brettsDesktop = 'C:\Users\baimi';
brettsMacbook = '/Users/bisraels';
clairesMacbook = '/Users/clairealbrecht';

% Name the computers and create string for  terminal string
computer_IDs = deblank(char(brettsDesktop, brettsMacbook, clairesMacbook));
computer_names = char('brettsDesktop', 'brettsMacbook', 'clairesMacbook');
% mode_str  = char('computer_baimi_mode', 'computer_bisraels_mode','computer_claire_mode');

numComputers = numel(computer_IDs(:,1));
n = 1;
while n <= numComputers
    %   TF = CONTAINS(STR,PATTERN) returns 1 (true) if STR contains PATTERN,
    %   and returns 0 (false) otherwise.
    terminalID = deblank(computer_IDs(n,:));
    if contains(wd, terminalID) == 1   % Find identifying piece of working directory
        user = deblank(computer_names(n,:));    % Associate the user name
        break;
    end
    if n == numComputers
    disp('Cannot identify computer. Setting user to ''general''.');
    user = 'general';           % If not one of the normal computers - general computer
    computer_general_mode = 1;
% % %     useFigurePosnMode = 0;
    end
end


% if verboseMode == 1 
%         disp(strcat('Turning on: ',computer_terminal_str));
%         gen_str = 'general';
%         if strcmp(user(1:5),gen_str(1:5)) == 1
% %             useFigurePosnMode = 0;
%             disp('Setting useFigurePosnMode to zero');
%         end


end
 
    
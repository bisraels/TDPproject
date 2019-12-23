function [computer_terminal_str, terminalID] = computerMode()
wd = pwd;

% verboseMode = 1;

% Define identifying strings in each computer's command line
baimi = 'C:\Users\baimi\';
bisraels = '/Users/bisraels/';
claire = '/Users/clairealbrecht/';

% Name the computers and create string for  terminal string
computer_IDs = char(baimi, bisraels, claire);
computer_names = char('baimi', 'bisraels', 'claire');
% mode_str  = char('computer_baimi_mode', 'computer_bisraels_mode','computer_claire_mode')

numComputers = numel(computer_IDs(:,1));
for i = 1:numComputers
    if contains(wd, deblank(computer_IDs(i,:)))   % Find identifying piece of working directory
        user = deblank(computer_names(i,:));  % Associate the user name
        terminalID = deblank(computer_IDs(i,:));
        disp(['Assigning terminalID as ' terminalID]);
    else
        user = 'general';           % If not one of the normal computers - general computer
        computer_general_mode = 1;
        useFigurePosnMode = 0;
    end
end

computer_terminal_str = ['computer_' user '_mode'];  % Define the output


% if verboseMode == 1 
%         disp(strcat('Turning on: ',computer_terminal_str));
%         gen_str = 'general';
%         if strcmp(user(1:5),gen_str(1:5)) == 1
% %             useFigurePosnMode = 0;
%             disp('Setting useFigurePosnMode to zero');
%         end
end
 
    
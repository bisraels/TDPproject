%__________________________________________________________________________
% AUTHOR:   Claire Albrecht
%
% CREATED:  November 2019
%
% PURPOSE:  Import solutions for weighting coefficients from mathematica
%           and correct format for our code.
%
% INPUT:    .csv file of mathematica solutions (wcoef_7state_us.nb)
%
% OUTPUT:   .mat file of symbolic cn_m in terms of vn_m
%
% MODIFICATION LOG:
%__________________________________________________________________________


tic

% Import .csv file with no heading
sol = readtable('wcoef_7state_varSol.csv','ReadVariableNames',0);
% First column tells you the variable name (for bookkeeping purposes only  - don't use this!)
% Second column is the solution we want.

% convert table into cell
sol_cell = table2cell(sol);

disp(['Replacing the "us" in string by "_" .Converting the long string into a sym 49 times.']);
% Replace the 'us' in the vector component variables with an UnderScore.
tic
c1_1 = sol_cell{1,2};
c1_1 = strrep(c1_1, 'us','_');
c1_1 = str2sym(c1_1);
elapsedTime = toc;
disp(['Time to do 1/49th = ' num2str(elapsedTime)]);
return

c1_2 = sol_cell{2,2};
c1_2 = strrep(c1_2, 'us','_');
c1_2 = str2sym(c1_2);

c1_3 = sol_cell{3,2};
c1_3 = strrep(c1_3, 'us','_');
c1_3 = str2sym(c1_3);

c1_4 = sol_cell{4,2};
c1_4 = strrep(c1_4, 'us','_');
c1_4 = str2sym(c1_4);

c1_5 = sol_cell{5,2};
c1_5 = strrep(c1_5, 'us','_');
c1_5 = str2sym(c1_5);

c1_6 = sol_cell{6,2};
c1_6 = strrep(c1_6, 'us','_');
c1_6 = str2sym(c1_6);

c1_7 = sol_cell{7,2};
c1_7 = strrep(c1_7, 'us','_');
c1_7 = str2sym(c1_7);

elapsedTime = toc;
disp(['Time to do 7/49th = ' num2str(elapsedTime)]);

c2_1 = sol_cell{8,2};
c2_1 = strrep(c2_1, 'us','_');
c2_1 = str2sym(c2_1);

c2_2 = sol_cell{9,2};
c2_2 = strrep(c2_2, 'us','_');
c2_2 = str2sym(c2_2);

c2_3 = sol_cell{10,2};
c2_3 = strrep(c2_3, 'us','_');
c2_3 = str2sym(c2_3);

c2_4 = sol_cell{11,2};
c2_4 = strrep(c2_4, 'us','_');
c2_4 = str2sym(c2_4);

c2_5 = sol_cell{12,2};
c2_5 = strrep(c2_5, 'us','_');
c2_5 = str2sym(c2_5);

c2_6 = sol_cell{13,2};
c2_6 = strrep(c2_6, 'us','_');
c2_6 = str2sym(c2_6);

c2_7 = sol_cell{14,2};
c2_7 = strrep(c2_7, 'us','_');
c2_7 = str2sym(c2_7);

elapsedTime = toc;
disp(['Time to do 2/7th = ' num2str(elapsedTime)]);

c3_1 = sol_cell{15,2};
c3_1 = strrep(c3_1, 'us','_');
c3_1 = str2sym(c3_1);

c3_2 = sol_cell{16,2};
c3_2 = strrep(c3_2, 'us','_');
c3_2 = str2sym(c3_2);

c3_3 = sol_cell{17,2};
c3_3 = strrep(c3_3, 'us','_');
c3_3 = str2sym(c3_3);

c3_4 = sol_cell{18,2};
c3_4 = strrep(c3_4, 'us','_');
c3_4 = str2sym(c3_4);

c3_5 = sol_cell{19,2};
c3_5 = strrep(c3_5, 'us','_');
c3_5 = str2sym(c3_5);

c3_6 = sol_cell{20,2};
c3_6 = strrep(c3_6, 'us','_');
c3_6 = str2sym(c3_6);

c3_7 = sol_cell{21,2};
c3_7 = strrep(c3_7, 'us','_');
c3_7 = str2sym(c3_7);

elapsedTime = toc;
disp(['Time to do 3/7th = ' num2str(elapsedTime)]);


c4_1 = sol_cell{22,2};
c4_1 = strrep(c4_1, 'us','_');
c4_1 = str2sym(c4_1);

c4_2 = sol_cell{23,2};
c4_2 = strrep(c4_2, 'us','_');
c4_2 = str2sym(c4_2);

c4_3 = sol_cell{24,2};
c4_3 = strrep(c4_3, 'us','_');
c4_3 = str2sym(c4_3);

c4_4 = sol_cell{25,2};
c4_4 = strrep(c4_4, 'us','_');
c4_4 = str2sym(c4_4);

c4_5 = sol_cell{26,2};
c4_5 = strrep(c4_5, 'us','_');
c4_5 = str2sym(c4_5);

c4_6 = sol_cell{27,2};
c4_6 = strrep(c4_6, 'us','_');
c4_6 = str2sym(c4_6);

c4_7 = sol_cell{28,2};
c4_7 = strrep(c4_7, 'us','_');
c4_7 = str2sym(c4_7);

elapsedTime = toc;
disp(['Time to do 4/7th = ' num2str(elapsedTime)]);



c5_1 = sol_cell{29,2};
c5_1 = strrep(c5_1, 'us','_');
c5_1 = str2sym(c5_1);

c5_2 = sol_cell{30,2};
c5_2 = strrep(c5_2, 'us','_');
c5_2 = str2sym(c5_2);

c5_3 = sol_cell{31,2};
c5_3 = strrep(c5_3, 'us','_');
c5_3 = str2sym(c5_3);

c5_4 = sol_cell{32,2};
c5_4 = strrep(c5_4, 'us','_');
c5_4 = str2sym(c5_4);

c5_5 = sol_cell{33,2};
c5_5 = strrep(c5_5, 'us','_');
c5_5 = str2sym(c5_5);

c5_6 = sol_cell{34,2};
c5_6 = strrep(c5_6, 'us','_');
c5_6 = str2sym(c5_6);

c5_7 = sol_cell{35,2};
c5_7 = strrep(c5_7, 'us','_');
c5_7 = str2sym(c5_7);

elapsedTime = toc;
disp(['Time to do 5/7th = ' num2str(elapsedTime)]);


c6_1 = sol_cell{36,2};
c6_1 = strrep(c6_1, 'us','_');
c6_1 = str2sym(c6_1);

c6_2 = sol_cell{37,2};
c6_2 = strrep(c6_2, 'us','_');
c6_2 = str2sym(c6_2);

c6_3 = sol_cell{38,2};
c6_3 = strrep(c6_3, 'us','_');
c6_3 = str2sym(c6_3);

c6_4 = sol_cell{39,2};
c6_4 = strrep(c6_4, 'us','_');
c6_4 = str2sym(c6_4);

c6_5 = sol_cell{40,2};
c6_5 = strrep(c6_5, 'us','_');
c6_5 = str2sym(c6_5);

c6_6 = sol_cell{41,2};
c6_6 = strrep(c6_6, 'us','_');
c6_6 = str2sym(c6_6);

c6_7 = sol_cell{42,2};
c6_7 = strrep(c6_7, 'us','_');
c6_7 = str2sym(c6_7);

elapsedTime = toc;
disp(['Time to do 6/7th = ' num2str(elapsedTime)]);


c7_1 = sol_cell{43,2};
c7_1 = strrep(c7_1, 'us','_');
c7_1 = str2sym(c7_1);

c7_2 = sol_cell{44,2};
c7_2 = strrep(c7_2, 'us','_');
c7_2 = str2sym(c7_2);

c7_3 = sol_cell{45,2};
c7_3 = strrep(c7_3, 'us','_');
c7_3 = str2sym(c7_3);

c7_4 = sol_cell{46,2};
c7_4 = strrep(c7_4, 'us','_');
c7_4 = str2sym(c7_4);

c7_5 = sol_cell{47,2};
c7_5 = strrep(c7_5, 'us','_');
c7_5 = str2sym(c7_5);

c7_6 = sol_cell{48,2};
c7_6 = strrep(c7_6, 'us','_');
c7_6 = str2sym(c7_6);

c7_7 = sol_cell{49,2};
c7_7 = strrep(c7_7, 'us','_');
c7_7 = str2sym(c7_7);

elapsedTime = toc;
disp(['Time to do 7/7th = ' num2str(elapsedTime)]);




% Save variables as a .mat file
save('wcoef_7state_cn_m.mat',...
    'c1_1','c1_2','c1_3','c1_4','c1_5','c1_6','c1_7',...
    'c2_1','c2_2','c2_3','c2_4','c2_5','c2_6','c2_7',...
    'c3_1','c3_2','c3_3','c3_4','c3_5','c3_6','c3_7',...
    'c4_1','c4_2','c4_3','c4_4','c4_5','c4_6','c4_7',...
    'c5_1','c5_2','c5_3','c5_4','c5_5','c5_6','c5_7',...
    'c6_1','c6_2','c6_3','c6_4','c6_5','c6_6','c6_7',...
    'c7_1','c7_2','c7_3','c7_4','c7_5','c7_6','c7_7')






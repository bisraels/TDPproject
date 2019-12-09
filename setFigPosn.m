function figPosn = setFigPosn(terminalID)

% Each  row of figPosn corresponds to the position of the figure for that
% row number.


baimi = 'C:\Users\baimi\';
bisraels = '/Users/bisraels/';
claire = '/Users/clairealbrecht/';

if strcmp(terminalID(1:10), baimi(1:10)) == 1           % Work Desktop
%     fig_1_posn = [1963 573 560 420];
%     fig_2_posn = [2535 571 560 420];
%     fig_3_posn = [3117 569 560 420];
%     fig_4_posn = [1963 48 560 420];
%     fig_5_posn = [3100 60 737 420];
%     fig_6_posn = [2534 48 560 420];
%     fig_7_posn = [1443 72 459 332];
    
    fig_1_posn = [1963; 573; 560; 420];
    fig_2_posn = [2535; 571; 560; 420];
    fig_3_posn = [3117; 569; 560; 420];
    fig_4_posn = [1963; 48; 560; 420];
    fig_5_posn = [3100; 60; 737; 420];
    fig_6_posn = [2534; 48; 560; 420];
    fig_7_posn = [1443; 72; 459; 332];
    figPosn = [fig_1_posn; fig_2_posn; fig_3_posn; fig_4_posn; fig_5_posn; fig_6_posn; fig_7_posn];
    
elseif strcmp(terminalID(1:10), bisraels(1:10)) == 1    % Macbook pro
%     fig_1_posn = [1 535 560 420];
%     fig_2_posn = [561 535 560 420];
%     fig_3_posn = [1121 535 560 420];
%     fig_4_posn = [1 41 560 420];
%     fig_5_posn = [562 41 842 420];
%     fig_6_posn = [1121 51 560 420];
%     fig_7_posn = [559 31 560 420];
    
    fig_1_posn = [1; 535; 560; 420];
    fig_2_posn = [561; 535; 560; 420];
    fig_3_posn = [1121; 535; 560; 420];
    fig_4_posn = [1; 41; 560; 420];
    fig_5_posn = [562; 41; 842; 420];
    fig_6_posn = [1121; 51; 560; 420];
    fig_7_posn = [559; 31; 560; 420];
    figPosn = [fig_1_posn; fig_2_posn; fig_3_posn; fig_4_posn; fig_5_posn; fig_6_posn; fig_7_posn];

elseif strcmp(terminalID(1:10), claire(1:10)) == 1
    %     fig_1_posn = [1 535 560 420];
    %     fig_2_posn = [561 535 560 420];
    %     fig_3_posn = [1121 535 560 420];
    %     fig_4_posn = [1 41 560 420];
    %     fig_5_posn = [562 41 842 420];
    %     fig_6_posn = [1121 51 560 420];
    %     fig_7_posn = [559 31 560 420];
    %     figPosn = [fig_1_posn; fig_2_posn; fig_3_posn; fig_4_posn; fig_5_posn; fig_6_posn; fig_7_posn];
    
    fig_1_posn = [1;438;508;367];
    fig_2_posn = [507;439;468;366];
    fig_3_posn = [976;438;465;367];
    fig_4_posn = [1;1;521;363];
    fig_5_posn = [1121; 51; 560; 420;];
    fig_6_posn = [915;1;526;363];
    fig_7_posn = [559; 31; 560; 420;];
    figPosn = [fig_1_posn, fig_2_posn, fig_3_posn, fig_4_posn, fig_5_posn, fig_6_posn, fig_7_posn];
    
end
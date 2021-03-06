% AUTHOR: Claire Albrecht
%
% CREATED:    July 2019
% MODIFIED:
%
% PURPOSE: Calculate probability in terms of rates (for TDP plots)
%
% Using the fact that at detailed balance (equilibrium) the rate in change
% of probability in each state is zero.
%               dP/dt = K * P = 0
% where K is the rate matrix and P is a vector of probabilities.
%
% The hard part is constructing the rate matrix, so this program will be
% written for states up to six. Since I cannot currently find a way to
% generalize the structure of the matrix, I will write them by hand here.
%
% MODIFICATIONS: August 2019 - included like to calculate evecs & evals 
%                                NOTE: Takes long time for general matrix,
%                                quicker once define some zeros.
%               August 2019  - added n = 4.1 to do the four state system
%               used in the Carey/Brett system with 0, 1, 2, 3 (labeling
%               (0, 1, 1', 2 states) and only 2 is connected to 3.
%               NOTE: This gives output of P#_Cn4 (with the C indicating
%               this is the carey system with reduced connections in the matrix)


% Number of states
n = 2;

syms k12 k13 k14 k15 k16 k21 k23 k24 k25 k26 k31 k32 k34 k35 k36 k41 k42 k43 k45 k46 k51 k52 k53 k54 k56 k61 k62 k63 k64 k65 
syms P1 P2 P3 P4 P5 P6

if n == 2
    K = [- k12 k21; k12 -k21];
 %   [V_n2, D_n2] = eig(K);
    P = sym('P',[1 2],'real');
    eqns = [mtimes(K,P') == 0; P1 + P2 == 1];
    
    solns = solve(eqns,[P1 P2]);
    P1_n2 = solns.P1
    P2_n2 = solns.P2

    save('sym_prob.mat','P1_n2','P2_n2')
    
end

if n == 3
    K = [(-k12 - k13) k21 k31; k12 (-k21 - k23) k32; k13 k23 (-k31 - k32)];
   % [V_n3, D_n3] = eig(K);
    P = sym('P',[1 3],'real');
    
    eqns = [mtimes(K,P') == 0; P1 + P2 + P3 == 1];
    
    solns = solve(eqns,[P1 P2 P3]);
    P1_n3 = solns.P1
    P2_n3 = solns.P2
    P3_n3 = solns.P3
    
    save('sym_prob.mat','P1_n3','P2_n3','P3_n3','-append')
end

if n == 4
    K = [(-k12 - k13 - k14) k21 k31 k41; k12 (-k21 - k23 - k24) k32 k42; k13 k23 (-k31 - k32 - k34) k43; k14 k24 k34 (-k41 - k42 - k43)];
  %  [V_n4, D_n4] = eig(K);
    P = sym('P',[1 4],'real');
    eqns = [mtimes(K,P') == 0; P1 + P2 + P3 + P4 == 1];
    
    solns = solve(eqns,[P1 P2 P3 P4]);
    P1_n4 = solns.P1
    P2_n4 = solns.P2
    P3_n4 = solns.P3
    P4_n4 = solns.P4
    
   % save('sym_prob.mat','P1_n4','P2_n4','P3_n4','P4_n4','V_n4','D_n4','-append')
   save('sym_prob.mat','P1_n4','P2_n4','P3_n4','P4_n4','V_n4','D_n4')
end


%--------------------------------------------------------------------------
% Use this for Carey model with states 0,1,2,3 (which label 0,1,1',2 states)
%--------------------------------------------------------------------------
if n == 4.1
    syms P0 k01 k02 k03 k10 k20 k30
% K = [-(k01 + k02 + 0), k10, k20, 0;...
%     k01, -(k10 + k12 + 0), k21, 0;...
%     k02, k12, -(k20 + k21 + k23), k32;...
%     0, 0, k23, -(0 + 0 + k32);];
K = [-(k01 + k02 + 0), k10, k20, 0;...
       k01, -(k10 + k12 + 0), k21, 0;...
       k02, k12, -(k20 + k21 + k23), k32;...
       0, 0, k23, -(0 + 0 + k32);];
%  [V_n4, D_n4] = eig(K);
    P = [P0; P1; P2; P3];
    assume(P0,'real')
    assume(P1,'real')
    assume(P2,'real')
    assume(P3,'real')
    
    % Detailed balance
    k02 = (k01 * k12 * k20)/(k10 * k21);
    
    % eqns = [mtimes(K,P) == 0; P0 + P1 + P2 + P3 == 1];
    eqns = [K * P == 0; P0 + P1 + P2 + P3 == 1];
    solns = solve(eqns,[P0 P1 P2 P3]);
    P0_Cn4 = solns.P0;
    P1_Cn4 = solns.P1;
    P2_Cn4 = solns.P1;
    P3_Cn4 = solns.P1;
    
    % save('sym_prob.mat','P0_Cn4','P1_Cn4','P2_Cn4','P3_Cn4','-append')
     save('sym_prob.mat','P0_Cn4','P1_Cn4','P2_Cn4','P3_Cn4')
end
%--------------------------------------------------------------------------



if n == 5
    K = [(-k12 - k13 - k14 - k15) k21 k31 k41 k51; k12 (-k21 - k23 - k24 - k25) k32 k42 k52; k13 k23 (-k31 - k32 - k34 - k35) k43 k53; k14 k24 k34 (-k41 - k42 - k43 - k45) k54; k15 k25 k35 k45 (-k51 - k52 - k53 - k54)];
  %  [V_n5, D_n5] = eig(K);
    P = sym('P',[1 5],'real');
    eqns = [mtimes(K,P') == 0; P1 + P2 + P3 + P4 + P5 == 1];
    
    solns = solve(eqns,[P1 P2 P3 P4 P5]);
    P1_n5 = solns.P1
    P2_n5 = solns.P2
    P3_n5 = solns.P3
    P4_n5 = solns.P4
    P5_n5 = solns.P5
        
    eigenvalue_n5 = eig(K);
    
    save('sym_prob.mat','P1_n5','P2_n5','P3_n5','P4_n5','P5_n5','V_n5','D_n5','-append')
end

% Note: n = 6 (and higher) is too much for my computer to do. Run it on computer down
% stairs.
if n == 6
    K = [(-k12 - k13 - k14 - k15 - k16) k21 k31 k41 k51 k61; k12 (-k21 - k23 - k24 - k25 - k26) k32 k42 k52 k62; k13 k23 (-k31 - k32 - k34 - k35 - k36) k43 k53 k63; k14 k24 k34 (-k41 - k42 - k43 - k45 -k46) k54 k64; k15 k25 k35 k45 (-k51 - k52 - k53 - k54 - k56) k65; k16 k26 k36 k46 k56 (-k61 - k62 - k63 - k64 - k65)];
 %   [V_n6, D_n6] = eig(K);
    P = sym('P',[1 6],'real');
    eqns = [mtimes(K,P') == 0; P1 + P2 + P3 + P4 + P5 + P6 == 1];
    
    solns = solve(eqns,[P1 P2 P3 P4 P5 P6]);
    P1_n6 = solns.P1
    P2_n6 = solns.P2
    P3_n6 = solns.P3
    P4_n6 = solns.P4
    P5_n6 = solns.P5
    P6_n6 = solns.P6
    
    save('sym_prob.mat','P1_n6','P2_n6','P3_n6','P4_n6','P5_n6','P6_n6','V_n6','D_n6','-append')

end

%%



%%

%x = load('temp.mat');
%y = load('sym_prob.mat');

%f1 = cell2mat(fieldnames(x));
%f2 = cell2mat(fieldnames(y));

% concatenate
%Ptot = [x.f1;y.f2]

%this combines the n6 results from downstairs with the n<6 calculations
%done here.



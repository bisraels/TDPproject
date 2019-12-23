function [histogram_FilePath, C2_FilePath, C4_FilePath] = loadFilePath_dataPlots(genAlgDataLocationPrefix, protein_str)

global verboseMode
%--------------------------------------------------------------------------
% HISTOGRAM:
%--------------------------------------------------------------------------
% fnamekeyword = '*_0p0uMgp32_001000us_GaussianFit_N3.mat';
% fNameKeyWord = '*_0p0uMgp32_001000us_FREThistogram.mat';
% fnamekeyword = '3p15mer_0p0uMgp32_000100us_FREThistogram.mat';
% fnamekeyword = '3p15mer_0p0uMgp32_000010us_FREThistogram.mat';
% fNameKeyWord = '3p15mer_0p5uMgp32_001000us_Neq038_FREThistogram.mat';

fNameKeyWord = '3p15mer_0p5uMgp32_001000us_FREThistogram.mat';

fileNames = dir([genAlgDataLocationPrefix filesep() fNameKeyWord]);
if isempty(fileNames) == 1
    error(['Could not Find any files with the name '  fNameKeyWord]);
end
histogram_FileName = fileNames(1).name;
histogram_FilePath = [genAlgDataLocationPrefix filesep() histogram_FileName];
if verboseMode ==  1
    if contains(histogram_FilePath(length(genAlgDataLocationPrefix):end),protein_str(1:4))&&contains(histogram_FilePath(length(genAlgDataLocationPrefix):end),protein_str(6:end))
        disp('We have the correct protein condition for histogram.')
    end
end

%--------------------------------------------------------------------------
% C2:
%--------------------------------------------------------------------------
% fNameKeyWord = '*_0p0uMgp32_tcfavg_2exp_fitresult.mat';
fNameKeyWord  = '3p15mer_0p5uMgp32_C2composite.mat';

fileNames = dir([genAlgDataLocationPrefix filesep() fNameKeyWord]);
if isempty(fileNames) == 1
    error(['Could not Find any files with the name '  fNameKeyWord]);
end
C2_FileName = fileNames(1).name;
C2_FilePath = [genAlgDataLocationPrefix filesep() C2_FileName];

if  verboseMode == 1
    if contains(C2_FilePath(length(genAlgDataLocationPrefix):end),protein_str(1:4))&&contains(C2_FilePath(length(genAlgDataLocationPrefix):end),protein_str(6:end))
        disp('We have the correct protein condition for C2.')
    end
end
%--------------------------------------------------------------------------
% C4:
%--------------------------------------------------------------------------
% fnamekeyword = '000010us_tau2-000000us_Neq035_fourptTCFavg.mat';

% fNameKeyWord = '*0p0uMgp32_tau2eq000000us_C4.mat';
fNameKeyWord = '3p15mer_0p5uMgp32_fourExpYoffset_Fit.mat';

fileNames = dir([genAlgDataLocationPrefix filesep() fNameKeyWord]);
if isempty(fileNames) == 1
    error(['Could not Find any files with the name '  fNameKeyWord]);
end
C4_FileName = fileNames(1).name;
C4_FilePath = [genAlgDataLocationPrefix filesep() C4_FileName];

if verboseMode == 1
    if contains(C4_FilePath(length(genAlgDataLocationPrefix):end),protein_str(1:4))&&contains(C4_FilePath(length(genAlgDataLocationPrefix):end),protein_str(6:end))
        disp('We have the correct protein condition for C4.')
    end
end


        
tic; % start timer
%% Clearing workspace and adding necessary files to path
clear; clc; close all;
restoredefaultpath;
addpath(genpath('lib/'));
rng('default');

codeVersion = 'v6.1';
codeDate = '2023_06_07';

%% This code performs 3 step CPCA on two sets of variables.
%
% Make sure to
% 1. have all your data in one .csv file.
% 2. encode all the missing data either using -99 or -88.
% 3. Have or install the 'Bioinformatics' and 'the Stat1istics and Machine
% learning' toolbxes.
%

%% Parameters you need to enter:
data = readcell('Data/current/IDs_Cybocs_current_no_yes_merged_CANTAB.csv'); % Insert your path/filename for your data here
    
Z_index_cell = {[14:31]}; 
G_index_cell = {[2:13]};
Sub_index_cell = {2:108}; % Enter the row numbers for subject groups on which CPCA analysis needs to be done as seen in .csv file

G_label = {'CANTAB'}; 
Z_label = {'CYBOCS'};
dataSetName = 'CYBOCS-in-Z__CANTAB-in-G_current__2-comps'; % Name of your dataset
% Function to rename labels
Sub_label = {'Subjects'};
reportsDir = sprintf('%s_Reports/', dataSetName);

%% Setting certain parameters for CPCA

params.nIter = 100; % No. times to split data in half
params.varimaxFlag = true;
params.p_val = 0.05;
params.BH_p_val = 0.05;
params.n_bootstrap = 100; % Shuffling values in Z and G to get reliability values for random data (to compare to real)
params.rValues = [0.1, 0.15, 0.2, 0.25,0.3,0.35,0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7];
params.addNumberToLabels = false;

% default settings - will not ask for inputs
params.default_flag = true;
params.default_struct.nComp = 2; % No. of comps
params.default_struct.sigType = 2; % Type of significance 1: p<=0.005, 2. p<=0.05 MC
params.default_struct.nDomCompLoadings = [6,6,6]; % Select the number of dominant component loadings for each component


%%  START OF SCRIPT

% Also change the Z, G, and Sub labels for each set that you enter above
% e.g. if your first set of Z variables are SSPI metrics, then change 'Z variable set 1' to 'SSPI'

%% Set missing values to NaN

data_num_only = data(2:end,1:end);
data_num_only(cellfun(@(s)ischar(s),data_num_only)) = {NaN};
data_num_only(cellfun(@(s)ismissing(s),data_num_only)) = {NaN};
data_num_only(cellfun(@(s)isequal(s,-99),data_num_only)) = {NaN};
data_num_only(cellfun(@(s)isequal(s,-88),data_num_only)) = {NaN};

data(2:end,1:end) = data_num_only;

%% Renames the variables
data(1,1:end) = strrep(data(1,1:end),'_',' ');


% Load in or create PLRP difference arrays
if isfile(strcat(reportsDir, "/rValues_cumulative.mat"))    % File exists.
     rValues_cumulative = load(strcat(reportsDir, "/rValues_cumulative.mat")).rValues_cumulative;
else                                                       % File does not exist.
    rValues_cumulative = [];
end
if isfile(strcat(reportsDir, "/PLRP_differences_CPCA.mat")) % File exists.
     PLRP_differences_CPCA = load(strcat(reportsDir, "/PLRP_differences_CPCA.mat")).PLRP_differences_CPCA;
else                                                       % File does not exist.
    PLRP_differences_CPCA = [];
end
if isfile(strcat(reportsDir, "/PLRP_differences_PCA.mat")) % File exists.
     PLRP_differences_PCA = load(strcat(reportsDir, "/PLRP_differences_PCA.mat")).PLRP_differences_PCA;
else                                                      % File does not exist.
    PLRP_differences_PCA = [];
end



% for r = params.rValues
for i = 1:length(params.rValues)

    r = params.rValues(i);
    params.r = r;

    %% Report settings

    import mlreportgen.dom.*;
    import mlreportgen.report.*;

    savePath = sprintf('%s_Reports/%s_Reports_r_%.2f/', dataSetName,dataSetName, r);
    mkdir(savePath);
    reportObjects = reportsInit(savePath, dataSetName, codeVersion, codeDate, r);

    %% Main interation

    for zi = 1:length(Z_index_cell)
        reportObjects = reportsInitChapter(reportObjects, zi, sprintf('Z matrix: %s',Z_label{zi}));

        for gi = 1:length(G_index_cell)
            reportObjects = reportsInitSectionI(reportObjects, [zi, gi], sprintf('G matrix: %s',G_label{gi}));

            for si = 1:length(Sub_index_cell)
                reportObjects = reportsInitSectionII(reportObjects, [zi, gi, si], '');
                params.page_header = sprintf('Z matrix: %s, G matrix: %s, Subject set: %s',Z_label{zi},G_label{gi},Sub_label{si});

                [reportObjects, CPCA_beh_struct(zi,gi,si), PCA_beh_struct(zi,gi,si), ...
                    CPCA_fullZ_struct(zi,gi,si), PCA_fullZ_struct(zi,gi,si)] = cpcaGenerateReport(reportObjects, [zi,gi,si], ...
                    data, Z_index_cell{zi}, G_index_cell{gi}, Sub_index_cell{si}, params);

                reportObjects = reportsAddSectionIITitles(reportObjects, [zi, gi, si], sprintf('Subject set: %s (n = %d)',...
                    Sub_label{si},size(CPCA_fullZ_struct(zi,gi,si).Z,1)));

                reportObjects = reportsAddSectionII(reportObjects, [zi, gi, si]);
                pause(1);
            end
            reportObjects = reportsAddSectionI(reportObjects, [zi, gi]);
        end
        reportObjects = reportsAddChapter(reportObjects, zi);
        close all;
    end

    save(fullfile(savePath,sprintf('%s_data.mat',dataSetName)),'CPCA_beh_struct','PCA_beh_struct','CPCA_fullZ_struct','PCA_fullZ_struct');

    inStruct.PCA_beh_struct = PCA_beh_struct;
    inStruct.CPCA_beh_struct = CPCA_beh_struct;
    inStruct.PCA_fullZ_struct = PCA_fullZ_struct;
    inStruct.CPCA_fullZ_struct = CPCA_fullZ_struct;
    inStruct.G_label = G_label;
    inStruct.Z_label = Z_label;
    inStruct.Sub_label = Sub_label;

    %% NEW: for PLRP p r plots
    rValues_cumulative(1,end+1) = r;
    PLRP_differences_CPCA(1,end+1) =  CPCA_beh_struct.PLRP_p0_001 - CPCA_beh_struct.PLRP_p0_05;
    PLRP_differences_PCA(1,end+1) = PCA_beh_struct.PLRP_p0_001 - PCA_beh_struct.PLRP_p0_05;
    

    [reportObjects] = reportsAddVarTables(reportObjects, inStruct);
    [reportObjects] = reportsAddSummary(reportObjects, inStruct);
    reportObjects = reportsClose(reportObjects);

    clear CPCA* reportObjects inStruct;
    
end


%% Sort PLRP arrays 
[rValues_cumulative, I] = sort(rValues_cumulative);
PLRP_differences_CPCA = PLRP_differences_CPCA(I);
PLRP_differences_PCA = PLRP_differences_PCA(I);

%% Plot PLRP difference over r values
plot(rValues_cumulative,PLRP_differences_CPCA);
hold on;
plot(rValues_cumulative, PLRP_differences_PCA);
title('PLRP(p=0.001) - PLRP(p=0.05) over r values');
legend({'CPCA','PCA'},'Location','southwest');
xlabel('r value');
ylabel('PLRP difference');
hold off;

saveas(gcf,strcat(reportsDir, '/PLRP_diff_over_r.png'));
save(strcat(reportsDir, "/rValues_cumulative.mat"), "rValues_cumulative");
save(strcat(reportsDir, "/PLRP_differences_CPCA.mat"), "PLRP_differences_CPCA");
save(strcat(reportsDir, "/PLRP_differences_PCA.mat"), "PLRP_differences_PCA");

 toc % end timer 

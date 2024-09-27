function [reportObjects, CPCA_struct_beh, PCA_struct_beh, CPCA_struct_fullZ, PCA_struct_fullZ] = cpcaGenerateReport(reportObjects, IdxNos, orig_data, Z_index, G_index, Sub_index, params)

import mlreportgen.dom.*;
import mlreportgen.report.*;

rptObj = reportObjects.SectionsII(IdxNos(1), IdxNos(2), IdxNos(3)).rptObj;

%% Get Z matrix

Z = cell2mat(orig_data(Sub_index,Z_index));
Zlabel = orig_data(1,Z_index);
if params.addNumberToLabels
    for ii = 1:length(Zlabel)
        Zlabel{ii} = sprintf('%d. %s',ii,Zlabel{ii});
    end
end

%% Get G Matrix

G = cell2mat(orig_data(Sub_index,G_index));
% G = [G ones(size(G,1),1)];
Glabel = orig_data(1,G_index);

%% Add numbers to the string label
if params.addNumberToLabels
    for ii = 1:length(Glabel)
        Glabel{ii} = sprintf('%d. %s',ii,Glabel{ii});
    end
end

%% Remove columns with zero variance and rows with NaNs

ZRejIdx = std(Z,[],1,"omitnan") == 0 | isnan(std(Z,[],1,"omitnan"));
Z(:,ZRejIdx) = [];
ZRejLabel = Zlabel(ZRejIdx);
Zlabel(ZRejIdx) = [];

GRejIdx = std(G,[],1,"omitnan") == 0 | isnan(std(G,[],1,"omitnan"));
G(:,GRejIdx) = [];
GRejLabel = Glabel(GRejIdx);
Glabel(GRejIdx) = [];

Z_orig = Z; G_orig = G;
[Z, G, SubRejIdx] = removeNaN(Z, G);
SubRejIds = find(SubRejIdx);

%% Find which Z and G variables are missing for each subject

subIds = orig_data(Sub_index,1);

fid = fopen('MissingSubjectData.txt','w');
for ii = SubRejIds'
    fprintf(fid,'Subject %d (ID = %d)\n',ii,subIds{ii});
    fprintf(fid,'Missing Z variables = ');
    zMissLabels = Zlabel(isnan(Z_orig(ii,:)));
    if length(zMissLabels)>1
        fprintf(fid,'%s, ',zMissLabels{1:end-1});
        fprintf(fid,'%s \n',zMissLabels{end});
    elseif ~isempty(zMissLabels)
        fprintf(fid,'%s \n',zMissLabels{end});
    else
        fprintf(fid,'\n');
    end

    fprintf(fid,'Missing G variables = ');
    gMissLabels = Glabel(isnan(G_orig(ii,:)));
    if length(gMissLabels)>1
        fprintf(fid,'%s, ',gMissLabels{1:end-1});
        fprintf(fid,'%s \n',gMissLabels{end});
    elseif ~isempty(gMissLabels)
        fprintf(fid,'%s \n',gMissLabels{end});
    else
        fprintf(fid,'\n');
    end

    fprintf(fid,'\n\n');
end
fclose(fid);

%%

Zdim = size(Z,2);
Gdim = size(G,2);
n_sub = size(Z,1);

if isempty(Z) || isempty(G)
    P = Text('This analysis is not possible, because there is no data in either Z or G.');
    rptObj = reportsAddObj(rptObj, P);
    return;
end

params.Zlabel = Zlabel;
params.Glabel = Glabel;
params.G_orig = G;
params.Z_orig = Z;
params.Zdim = Zdim;
params.Gdim = Gdim;
params.n_sub = n_sub;

%% Initial page for behavioral and fullZ CPCA report

H = Heading4(Text('Behavioral CPCA analysis')); 
rptObj = reportsAddObj(rptObj, H);

T = Text(sprintf('n = %d \n', n_sub)); 
rptObj = reportsAddObj(rptObj, T);

if isempty(SubRejIds)
    T = Text(sprintf('No subjects rejected due to NaN in a variable'));
else
    T = Text([sprintf('Subjects rejected due to NaN in a variable = '), sprintf('%d, ',SubRejIds(1:end-1)),sprintf('%d',SubRejIds(end))]);
end
rptObj = reportsAddObj(rptObj, T);

T = Text([sprintf('Variables in Z rejected due to zero variance = '), sprintf('%s,',ZRejLabel{:})]); 
rptObj = reportsAddObj(rptObj, T);
T = Text([sprintf('Variables in G rejected due to zero variance = '), sprintf('%s,',GRejLabel{:})]); 
rptObj = reportsAddObj(rptObj, T);

T = Text(sprintf('%s \n', params.page_header)); 
rptObj = reportsAddObj(rptObj, T);

T = Text(sprintf('Z-dim = %d \n', Zdim)); 
rptObj = reportsAddObj(rptObj, T);
T = Text([sprintf('List of Z variables = '), sprintf('%s, ',Zlabel{:})]); 
rptObj = reportsAddObj(rptObj, T);


T = Text(sprintf('G-dim = %d \n', Gdim));  
rptObj = reportsAddObj(rptObj, T);
T = Text([sprintf('List of G variables = '), sprintf('%s, ',Glabel{:})]); 
rptObj = reportsAddObj(rptObj, T);

br = PageBreak();
rptObj = reportsAddObj(rptObj, br);

%% Basic full Z CPCA

rng('default');
[rptObj(1,:), CPCA_struct_fullZ, PCA_struct_fullZ] = runFullZCpca(rptObj(1,:), params);

%% Behavioral CPCA - with iterations

rng('default');
[rptObj(2,:), CPCA_struct_beh, PCA_struct_beh] = runBehCpca(rptObj(2,:), params);

reportObjects.SectionsII(IdxNos(1), IdxNos(2), IdxNos(3)).rptObj = rptObj;


end

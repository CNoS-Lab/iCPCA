function [rpt_obj, CPCA_struct, PCA_struct] = runBehCpca(rpt_obj, params)

%% runBehCpca 
% 
% Function to run the formally-called Behavioral-CPCA, now called
% Iterative-CPCA
% 
% Note: To change formatting, lines where fontSizes are set are commented
% with 'font'
% 
% Script parameters

%% Imports
import mlreportgen.dom.*;
import mlreportgen.report.*;

Z = params.Z_orig;
G = params.G_orig;

Z = (Z - mean(Z,1,"omitnan"))./std(Z,[],1,"omitnan"); % Z-score the data

%% Perform CPCA

if params.default_flag
    nComp = params.default_struct.nComp;
else
    nComp = [];
end

params_cpca = params;
params_cpca.cpca_pca_flag = 'cpca'; % Run CPCA

params_pca = params;
params_pca.cpca_pca_flag = 'pca'; % Run PCA

CPCA_struct = runIterCPCA(Z, G, nComp, params_cpca);
nComp = CPCA_struct.nComp;

permutePCorrRelPCA = NaN([size(CPCA_struct.mean_PCorr),params.n_bootstrap]);
permutePCorrRelCPCA = NaN([size(CPCA_struct.mean_PCorr),params.n_bootstrap]);

for i = 1:params.n_bootstrap
    G_temp = G(randperm(size(G,1)),:); Z_temp = Z;

    temp_struct = runIterCPCA(Z_temp, G_temp, nComp, params_cpca);
    permutePCorrRelCPCA(:,:,i) = temp_struct.mean_PCorr_reliable;

    temp_struct = runIterCPCA(Z_temp, G_temp, nComp, params_pca);
    permutePCorrRelPCA(:,:,i) = temp_struct.mean_PCorr_reliable;

    clc; disp(i);
end

CPCA_struct.predrel_cutoff_rel = quantile(abs(permutePCorrRelCPCA(:)), 1-params.p_val);
CPCA_struct.RelPredLoadMask = abs(CPCA_struct.mean_PCorr_reliable) >=  CPCA_struct.predrel_cutoff_rel;

CPCA_struct.nullDist = permutePCorrRelCPCA;
CPCA_struct.params = params_cpca;

% estimate the PLRP values for p=0.001 and p=0.05 for selection of r

CPCA_struct.PLRP_p0_001 = quantile(abs(permutePCorrRelCPCA(:)), 1-0.001);
CPCA_struct.PLRP_p0_05 = quantile(abs(permutePCorrRelCPCA(:)), 1-0.05);

%% Perform PCA

PCA_struct = runIterCPCA(Z, G, nComp, params_pca);

PCA_struct.predrel_cutoff_rel = quantile(abs(permutePCorrRelPCA(:)), 1-params.p_val);
PCA_struct.RelPredLoadMask = abs(PCA_struct.mean_PCorr_reliable) >=  PCA_struct.predrel_cutoff_rel;

PCA_struct.nullDist = permutePCorrRelPCA;
PCA_struct.params = params_pca;

% estimate the PLRP values for p=0.001 and p=0.05 for selection of r

PCA_struct.PLRP_p0_001 = quantile(abs(permutePCorrRelPCA(:)), 1-0.001);
PCA_struct.PLRP_p0_05 = quantile(abs(permutePCorrRelPCA(:)), 1-0.05);

%% Plot num of reliable predictor loadings - BH multiple comparisson correction

for ii = 1:size(CPCA_struct.PCorr,1)
    for jj = 1:size(CPCA_struct.PCorr,2)
        PCA_struct.pVales(ii,jj) = mean(abs(permutePCorrRelPCA(:)) > abs(PCA_struct.mean_PCorr_reliable(ii,jj)),"all");
        CPCA_struct.pVales(ii,jj) = mean(abs(permutePCorrRelCPCA(:)) > abs(CPCA_struct.mean_PCorr_reliable(ii,jj)),"all");
    end
end

[PCA_struct.pVales_corrected, PCA_struct.alpha_corrected] = multicmp(PCA_struct.pVales(:), 'fdr', 0.05);
PCA_struct.pVales_corrected = reshape(PCA_struct.pVales_corrected, size(PCA_struct.pVales));

[CPCA_struct.pVales_corrected, CPCA_struct.alpha_corrected] = multicmp(CPCA_struct.pVales(:), 'fdr', 0.05);
CPCA_struct.pVales_corrected = reshape(CPCA_struct.pVales_corrected, size(CPCA_struct.pVales));

PCA_struct.RelPredLoadMask_corrected = PCA_struct.pVales_corrected <= 0.05;
CPCA_struct.RelPredLoadMask_corrected = CPCA_struct.pVales_corrected <= 0.05;

%% add PCA plots to report

f = Figure;

rpt_obj_sub_pca = Section;
rpt_obj_sub_pca.Title = 'PCA analysis';
titleStr = 'PCA analysis';

PCA_struct.Res = SelectCompLoadingsBig(PCA_struct, params_pca);
[rpt_obj_sub_pca, f, PCA_struct] = addPlots(rpt_obj_sub_pca, f, PCA_struct, params_pca, titleStr);
add(rpt_obj(1), rpt_obj_sub_pca);

% Get summary
PCA_struct = addSummary(PCA_struct, params_pca);

%% add CPCA plots to report

rpt_obj_sub_cpca = Section;
rpt_obj_sub_cpca.Title = 'CPCA analysis';
titleStr = 'CPCA analysis';

CPCA_struct.Res = SelectCompLoadingsBig(CPCA_struct, params_cpca);
[rpt_obj_sub_cpca, f, CPCA_struct] = addPlots(rpt_obj_sub_cpca, f, CPCA_struct, params_cpca, titleStr);
add(rpt_obj(2), rpt_obj_sub_cpca);

% Get summary
CPCA_struct = addSummary(CPCA_struct, params_cpca);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function addPlotslideshow?
% 
% Adds all of the Scree plots, component loadings plots, predictor loadings
% plot, and dominant component loadings plot.
%
%% Params:
% rpt_obj:
% f: 
% cpca_struct
% params:
% titleStr: 
%
%% Returns: 
% rpt_obj: 
function [rpt_obj, f, cpca_struct] = addPlots(rpt_obj, f, cpca_struct, params, titleStr)

% Scree plot and variance tables
[rpt_obj, f] = addScreePlot(rpt_obj, f, cpca_struct, params, titleStr);
[rpt_obj, f] = addPercentVarTable(rpt_obj, f, cpca_struct, params, titleStr);
[rpt_obj, f] = addPercentVarSpiderPlot(rpt_obj, f, cpca_struct, params, titleStr);

% Comp Loadings
[rpt_obj, f] = addCompLoadingsPlot(rpt_obj, f, cpca_struct, params, titleStr);
[rpt_obj, f] = addCompLoadingsSpiderPlot(rpt_obj, f, cpca_struct, params, titleStr);

% Null Distribution for predictor loadings
[rpt_obj, f] = addPredLoadThreshPlot(rpt_obj, f, cpca_struct, params, titleStr);

% Predictor loadings
[rpt_obj, f] = addPredLoadingsPlot_wMultCompCorrection(rpt_obj, f, cpca_struct, params, titleStr);

% Dominant comp loadings
[rpt_obj, f, cpca_struct] = addDominantCompLoadPlot(rpt_obj, f, cpca_struct, params, titleStr);

end



%% function addScreePlot
% 
% generates the Scree plot which shows the ... TODO: Finish
%
%% Params: 
% rpt_obj:
% f:
% cpca_struct:
% params:
% titleStr:
% 
%% Returns: 
% rpt_obj:
% f:
function [rpt_obj, f] = addScreePlot(rpt_obj, f, cpca_struct, params, titleStr) %% Plot scree plot

screePlotFontSize = 6;

import mlreportgen.dom.*;
import mlreportgen.report.*;

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text(sprintf('%s - Scree Plot',titleStr))); add(rpt_obj, H);

fig = figure('Position',[0 0 1280 720],'visible','off'); % Scree plot
D = cpca_struct.DGH;
D = diag(D).^2;
D = D./sum(D)*100;
k = length(D);
plot(1:k,D(1:k),'-O');
xlabel('No. of Components');
ylabel('Percentage of variance');
title(sprintf('%s - Scree Plot',titleStr));
set(findall(gcf,'-property','FontSize'),'FontSize',screePlotFontSize); 

f(1) = Figure(fig); add(rpt_obj, f(1));

T = Text(sprintf('%d number of components were selected. \n\n',cpca_struct.nComp));
add(rpt_obj, T);

br = PageBreak();
add(rpt_obj, br);

end




%% function addPercentVarTable
%
% Params: 
% 
% Returns
function [rpt_obj, f] = addPercentVarTable(rpt_obj, f, cpca_struct, params, titleStr) %% Plot predictor loadings

import mlreportgen.dom.*;
import mlreportgen.report.*;

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text(sprintf('%s - Percent Variance explained for each Z variable',titleStr))); add(rpt_obj, H);

percent_var = mean(cat(2,cpca_struct.varianceTable.perVarEachS),2,"omitnan");
percent_var = percent_var./100;

T = cell2table([cpca_struct.params.Zlabel',sprintfc('%.2f',percent_var)]);
T.Properties.VariableNames = {'Zlabel', '% variance explained by G'};

Tbl = MATLABTable(T);
Tbl.TableEntriesStyle = {Height('20pt')};
add(rpt_obj, Tbl);
br = PageBreak();
add(rpt_obj, br);

end



%% function addPercentVarSpiderPlot
%
% Params: 
% 
% Returns
function [rpt_obj, f] = addPercentVarSpiderPlot(rpt_obj, f, cpca_struct, params, titleStr)

spiderPlotFontSize = 6;

import mlreportgen.dom.*;
import mlreportgen.report.*;

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text(sprintf('%s - Percent Variance explained for each Z variable - Spider Plot',titleStr))); add(rpt_obj, H);

% Get Z labels
Z_label_str = cpca_struct.params.Zlabel;

fig = figure('Position',[100 100 600 600],'visible','off','color','w');

% plot percent variance explained
idx = reshape(1:64,[8,8])';
idx = idx(2:6,2:6);
subplot(8,8,idx(:));

percent_var = mean(cat(2,cpca_struct.varianceTable.perVarEachS),2,"omitnan");
percent_var = percent_var./100;

nZ = size(cpca_struct.Z,2);
for ii = 1:length(Z_label_str)
    Z_label_str{ii} = [Z_label_str{ii}, ' ', sprintf('%.2f',percent_var(ii))];
end

% Plot the spider plot
if strcmp(params.cpca_pca_flag,'pca')
    spider_plot(percent_var', 'AxesLabels', Z_label_str, ...
        'AxesLimits', repmat([0,1],nZ,1)','Direction', 'roseplot',...
        'Color', [1 0 0], 'BackgroundColor' , 'w', 'AxesDisplay', 'one',...
        'FillOption', {'on'}, 'FillTransparency', 0.5, 'AxesOffset', 0); % Spider plot
elseif strcmp(params.cpca_pca_flag,'cpca')
    spider_plot([ones(size(percent_var)),percent_var,percent_var]', 'AxesLabels', Z_label_str, ...
        'AxesLimits', repmat([0,1.001],nZ,1)','Direction', 'roseplot',...
        'Color', [1 0 0; 0 1 1; 0 1 0], 'BackgroundColor' , 'w', 'AxesDisplay', 'one',...
        'FillOption', {'on'}, 'FillTransparency', 0.5, 'AxesOffset', 0); % Spider plot
else
    error('This should not happen! params.cpca_pca_flag should be either "pca" or "cpca"');
end

title('Percent vaiance explained');
set(findall(gcf,'-property','FontSize'),'FontSize', spiderPlotFontSize); 

f(end+1) = Figure(fig);
f(end).Scaling = 'none';
add(rpt_obj, f(end));
br = PageBreak();
add(rpt_obj, br);

end



%% function addCompLoadingsPlot
%
%& Params: 
%
% Returns:
function [rpt_obj, f] = addCompLoadingsPlot(rpt_obj, f, cpca_struct, params, titleStr) %% Plot component loadings

compLoadingsFontSize = 8;

import mlreportgen.dom.*;
import mlreportgen.report.*;

for ii = 1:cpca_struct.nComp
    xStr{ii} = sprintf('C%d',ii);
end

H = Heading4(Text(sprintf('%s - Component Loadings',titleStr)));
add(rpt_obj, H);

nz = 10;
temp = cpca_struct.loadings_VD_sqrtN_GH;
n_plots = ceil(size(temp,1)/nz);

for ii = 1:n_plots
    T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
    fig = figure('visible','off');
    if ii == n_plots
        h = plotLoadings(temp(nz*(ii-1)+1:end,:),xStr,params.Zlabel(nz*(ii-1)+1:end));
    else
        h = plotLoadings(temp(nz*(ii-1)+1:nz*ii,:),xStr,params.Zlabel(nz*(ii-1)+1:nz*ii));
    end
    h.Colormap = flipud(cbrewer('div','RdBu',64,'spline'));
    h.ColorLimits = [-1*max(abs(temp(:))) max(abs(temp(:)))];
    if n_plots == 1
        title(sprintf('%s - Component Loadings',titleStr));
    else
        title(sprintf('%s - Component Loadings plot no. %d of %d',titleStr,ii,n_plots));
    end
    %set(fig, 'Position',[0 0 500 400]);
    set(fig, 'Position',[0 0 700 600]);
    set(findall(gcf,'-property','FontSize'),'FontSize', compLoadingsFontSize); 
    f(end+1) = Figure(fig);
    f(end).Scaling = 'none';
    add(rpt_obj, f(end));
    br = PageBreak();
    add(rpt_obj, br);
end

end



%% function addCompLoadingsSpiderPlot
%
% Params: 
% 
% Returns
function [rpt_obj, f] = addCompLoadingsSpiderPlot(rpt_obj, f, cpca_struct, params, titleStr)

import mlreportgen.dom.*;
import mlreportgen.report.*;

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text(sprintf('%s - Component Loadings - Spider Plot',titleStr))); add(rpt_obj, H);

% Get Z labels
Z_label_str = cpca_struct.params.Zlabel;

fig = figure('Position',[100 100 700 700],'visible','off','color','w');

% plot percent variance explained
idx = reshape(1:64,[8,8])';
idx = idx(2:6,2:6);
subplot(8,8,idx(:));

data = cpca_struct.loadings_VD_sqrtN_GH;
ax_lim = 1;
data = [ax_lim*ones(size(data,1),1),data];

nZ = size(cpca_struct.Z,2);
nComp = cpca_struct.nComp;

% Plot the spider plot
if strcmp(params.cpca_pca_flag,'pca')
    color_mat = [1 0 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 0 0 0; 1 1 1];
elseif strcmp(params.cpca_pca_flag,'cpca')
    color_mat = [0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 0 0 0; 1 1 1];
else
    error('This should not happen! params.cpca_pca_flag should be either "pca" or "cpca"');
end

spider_plot(abs(data)', 'AxesLabels', Z_label_str, 'AxesLimits', repmat([0,ax_lim],nZ,1)',...
    'Direction', 'roseplot','AxesOffset', 0, 'AxesDisplay', 'one', 'Color', color_mat(1:nComp+1,:), ...
    'FillOption', 'on', 'FillTransparency', 0.5); % Spider plot

title(sprintf('%s: component loadings',titleStr));
legend([{''},sprintfc('C%d',1:nComp)], 'Location', 'southoutside');

set(findall(gcf,'-property','FontSize'),'FontSize',8); 
f(end+1) = Figure(fig);
f(end).Scaling = 'none';
add(rpt_obj, f(end));
br = PageBreak();
add(rpt_obj, br);

end



%% function addPredLoadThreshPlot
%
% Params: 
% 
% Returns
function [rpt_obj, f] = addPredLoadThreshPlot(rpt_obj, f, cpca_struct, params, titleStr) %% Predictor reliability threshold

predLoadThreshFontSize = 12;

import mlreportgen.dom.*;
import mlreportgen.report.*;

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text(sprintf('%s - Predictor loading reliability proportions threshold - p <= 0.005', titleStr)));
add(rpt_obj, H);

predRel = abs(cpca_struct.nullDist(:));
fig = figure('visible','off');
histogram(predRel);
hold on;
clear g;
g(1) = vline(quantile(predRel(:), 1-0.05)); g(1).LineWidth = 2; g(1).Color = 'r';
g(2) = vline(quantile(predRel(:), 1-0.01)); g(2).LineWidth = 2; g(2).Color = 'm';
g(3) = vline(quantile(predRel(:), 1-0.005)); g(3).LineWidth = 2; g(3).Color = 'g';
g(4) = vline(quantile(predRel(:), 1-0.001)); g(4).LineWidth = 2; g(4).Color = 'b';
legend(g,{'p = 0.05','p = 0.01','p = 0.005','p = 0.001'});
title(sprintf('%s - Predictor loading threshold', titleStr));
set(fig, 'Position',[0 0 600 500]);
set(findall(gcf,'-property','FontSize'),'FontSize', predLoadThreshFontSize);

T = Text(sprintf('Estimated Predictor loading threshold = %4f (at p <= %.4f) \n', cpca_struct.predrel_cutoff_rel,params.p_val));
add(rpt_obj, T);

f(end+1) = Figure(fig);
f(end).Scaling = 'none';
add(rpt_obj, f(end));

br = PageBreak();
add(rpt_obj, br);

end



%% function addPredLoadingsPlot_wMultCompCorrection
%
% generates the predictor loadings with the 'multiple Comparisons issue
% corrected for. 
%
% Params: 
% 
% Returns
function [rpt_obj, f] = addPredLoadingsPlot_wMultCompCorrection(rpt_obj, f, cpca_struct, params, titleStr) %% Plot predictor loadings

PredLoadThresh_MultCompCorrected_FontSize = 8;

import mlreportgen.dom.*;
import mlreportgen.report.*;

for ii = 1:cpca_struct.nComp
    xStr{ii} = sprintf('C%d',ii);
end

H = Heading4(Text(sprintf('%s - Predictor Loadings and Predictor Loading reliability proportions (Multiple comparison corrected)',titleStr)));
add(rpt_obj, H);

nz = 10;
temp = cpca_struct.mean_PCorr;
temp_plrp = cpca_struct.mean_PCorr_reliable;
n_plots = ceil(size(temp,1)/nz);

for ii = 1:n_plots
    T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
    fig = figure('visible','off');
    if ii == n_plots
        h = plotLoadingsPLRP4(temp(nz*(ii-1)+1:end,:),temp_plrp(nz*(ii-1)+1:end,:),xStr,params.Glabel(nz*(ii-1)+1:end),cpca_struct.pVales_corrected(nz*(ii-1)+1:end,:), [0.05,0.01,0.001], cpca_struct.predrel_cutoff_rel);
    else
        h = plotLoadingsPLRP4(temp(nz*(ii-1)+1:nz*ii,:),temp_plrp(nz*(ii-1)+1:nz*ii,:),xStr,params.Glabel(nz*(ii-1)+1:nz*ii),cpca_struct.pVales_corrected(nz*(ii-1)+1:nz*ii,:), [0.05,0.01,0.001], cpca_struct.predrel_cutoff_rel);
    end
    h.Colormap = flipud(cbrewer('div','RdBu',64,'spline'));
    h.CLim = [-0.3 0.3];
    if n_plots == 1
        title(sprintf('%s - Predictor Loadings (Predictor Loading Reliability Proportion, Multiple comparison corrected p-value)',titleStr));
    else
        title(sprintf('%s - Predictor Loadings (Predictor Loading Reliability Proportion, Multiple comparison corrected p-value) plot no. %d of %d ',titleStr,ii,n_plots));
    end
    set(fig, 'Position',[0 0 600 600]);
    set(findall(gcf,'-property','FontSize'),'FontSize',PredLoadThresh_MultCompCorrected_FontSize);
    f(end+1) = Figure(fig);
    f(end).Scaling = 'none';
    add(rpt_obj, f(end));
    T = Text(sprintf(['* next to reliability values = significance based on uncorrected p<=%.4f \n ...' ...
        '*,**,*** next to p-values = significance based on corrected p<=0.05,0.01,0.001 \n\n'],params.p_val));
    add(rpt_obj, T);
    br = PageBreak();
    add(rpt_obj, br);
end

end




%% function addDominantCompPlot
%
% Params: 
% 
% Returns
function [rpt_obj, f, cpca_struct] = addDominantCompLoadPlot(rpt_obj, f, cpca_struct, params, titleStr) %% Plot predictor loadings

DomCompLoadingsFontSize = 6;

import mlreportgen.dom.*;
import mlreportgen.report.*;

H = Heading4(Text(sprintf('%s - Dominant Component Loadings',titleStr)));
add(rpt_obj, H);

comp_rel = find(any(cpca_struct.RelPredLoadMask | cpca_struct.RelPredLoadMask_corrected,1));

meanPL_rel = cat(1,cpca_struct.Res.meanPCorr_rel);
meanPL_rel_corrected = cat(1,cpca_struct.Res.meanPCorr_rel_corrected);

[meanPL_rel_sort, meanPL_rel_sortIdx] = sort(meanPL_rel,1,'ascend');
[meanPL_rel_corrected_sort, meanPL_rel_corrected_sortIdx] = sort(meanPL_rel_corrected,1,'ascend');

if ~isempty(comp_rel)
    fig = figure('visible','on');
    subplot(1,2,1);
    temp = sum(abs(cpca_struct.mean_PCorr_reliable.*cpca_struct.RelPredLoadMask),1)./sum(cpca_struct.RelPredLoadMask,1);
    plot([meanPL_rel_sort; temp],'*-','LineWidth',2);
    legend(sprintfc('C%d',1:cpca_struct.nComp),'Location','bestoutside');
    title('Reliability based on p<=0.005');
    ylabel('Mean of reliable Pred. Load. Rel. Prop.');
    xlabel('Variable dominance rank #');

    subplot(1,2,2);
    temp = sum(abs(cpca_struct.mean_PCorr_reliable.*cpca_struct.RelPredLoadMask_corrected),1)./sum(cpca_struct.RelPredLoadMask_corrected,1);
    plot([meanPL_rel_corrected_sort; temp],'*-','LineWidth',2);
    legend(sprintfc('C%d',1:cpca_struct.nComp),'Location','bestoutside');
    title('Reliability based on p<=0.05, corrected for multiple comparison');
    ylabel('Mean of reliable Pred. Load. Rel. Prop.');
    xlabel('Variable dominance rank #');

    sgtitle(sprintf('%s - Dominant Component Loadings',titleStr));
    set(fig, 'Position',[0 0 1280 720]);
    set(findall(gcf,'-property','FontSize'),'FontSize',15);
end

for comp = comp_rel
    T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
    fig = figure('visible','off');

    subplot(1,2,1);
    temp = sum(abs(cpca_struct.mean_PCorr_reliable(:,comp).*cpca_struct.RelPredLoadMask(:,comp)),1)./sum(cpca_struct.RelPredLoadMask(:,comp),1);
    plot([meanPL_rel_sort(:,comp); temp],'*-','LineWidth',2);
    title(sprintf('Comp %d - Reliability based on p<=0.005',comp));
    ylabel('Mean of reliable Pred. Load. Rel. Prop.');
    xlabel('Variable dominance rank #');
    xticklabels([sprintfc('%d',meanPL_rel_sortIdx(:,comp)'),'all']);

    subplot(1,2,2);
    temp = sum(abs(cpca_struct.mean_PCorr_reliable(:,comp).*cpca_struct.RelPredLoadMask_corrected(:,comp)),1)./sum(cpca_struct.RelPredLoadMask_corrected(:,comp),1);
    plot([meanPL_rel_corrected_sort(:,comp); temp],'*-','LineWidth',2);
    title(sprintf('Comp %d - Reliability based on p<=0.05, corrected for multiple comparison',comp));
    ylabel('Mean of reliable Pred. Load. Rel. Prop.');
    xlabel('Variable dominance rank #');
    xticklabels([sprintfc('%d',meanPL_rel_corrected_sortIdx(:,comp)'),'all']);

    set(fig, 'Position',[100 100 500 300]);
    set(findall(gcf,'-property','FontSize'),'FontSize',DomCompLoadingsFontSize);

    f(end+1) = Figure(fig);
    f(end).Scaling = 'none';
    add(rpt_obj, f(end));

    T = Text([sprintf('Xlabel for 1st plot, Mean Pred. Load. for p<=0.005,comp %d = ',comp),...
        sprintf('%d, ', meanPL_rel_sortIdx(1:end-1,comp)),sprintf('%d',meanPL_rel_sortIdx(end,comp))]);
    add(rpt_obj, T);

    T = Text([sprintf('Xlabel for 2nd plot, Mean Pred. Load. for corrected p<=0.05 ,comp %d = ',comp),...
        sprintf('%d, ', meanPL_rel_corrected_sortIdx(1:end-1,comp)),sprintf('%d',meanPL_rel_corrected_sortIdx(end,comp))]);
    add(rpt_obj, T);

    br = PageBreak();
    add(rpt_obj, br);
end

cpca_struct.DomCompLoad.meanPL_rel = meanPL_rel;
cpca_struct.DomCompLoad.meanPL_rel_corrected = meanPL_rel_corrected;
cpca_struct.DomCompLoad.meanPL_rel_sort = meanPL_rel_sort;
cpca_struct.DomCompLoad.meanPL_rel_sortIdx = meanPL_rel_sortIdx;
cpca_struct.DomCompLoad.meanPL_rel_corrected_sor = meanPL_rel_corrected_sort;
cpca_struct.DomCompLoad.meanPL_rel_corrected_sortIdx = meanPL_rel_corrected_sortIdx;

end



%% function addSummary
%
% Params: 
% 
% Returns
function [cpca_struct] = addSummary(cpca_struct, params)

Add_more_summary_flag = true;
summ_no = 1;
cpca_struct.Summary = struct();
summTypeStr = {'LEFT plot (p<=0.005)','RIGHT plot (p<=0.05, MC)'};

comp_rel = find(any(cpca_struct.RelPredLoadMask | cpca_struct.RelPredLoadMask_corrected,1));

while Add_more_summary_flag && ~isempty(comp_rel)

    fprintf('Summary %d: \n',summ_no);

    if params.default_flag
        summaryType = params.default_struct.sigType;
    else
        summaryType = input(sprintf('Iterative %s analysis: Enter the type of significance you would like to use - \n1 (for p<=0.005, LEFT PLOT) or \n2 (for p<=0.05 corrected for multiple comparisons, RIGHT PLOT) \nEnter here :',upper(params.cpca_pca_flag)));
    end

    if ~(isscalar(summaryType) && (summaryType == 1 || summaryType == 2))
        disp('Please enter type as either 1 or 2.');
        continue;
    end

    if summaryType == 1
        comp_rel_summ = find(any(cpca_struct.RelPredLoadMask,1));
    else
        comp_rel_summ = find(any(cpca_struct.RelPredLoadMask_corrected,1));
    end

    for c = 1:length(comp_rel_summ)
        cpca_struct.Summary(summ_no,c).CompNo = comp_rel_summ(c);
        cpca_struct.Summary(summ_no,c).SummaryType = summaryType;
        
        if params.default_flag
            opt_comp_no = params.default_struct.nDomCompLoadings(comp_rel_summ(c));
        else
            opt_comp_no = input(sprintf('Enter the number of dominant component loadings for comp %d (based on %s): ',comp_rel_summ(c),summTypeStr{summaryType}));
        end

        if summaryType == 1
            opt_comp_load_idx = cpca_struct.DomCompLoad.meanPL_rel_sortIdx(1:opt_comp_no,comp_rel_summ(c));
        elseif summaryType == 2
            opt_comp_load_idx = cpca_struct.DomCompLoad.meanPL_rel_corrected_sortIdx(1:opt_comp_no,comp_rel_summ(c));
        else
            error('This should never happen!');
        end

        for ii = 1:length(opt_comp_load_idx)
            cpca_struct.Summary(summ_no,c).CompLoadings(ii).Name = params.Zlabel{opt_comp_load_idx(ii)};
            cpca_struct.Summary(summ_no,c).CompLoadings(ii).Value = cpca_struct.loadings_VD_sqrtN_GH(opt_comp_load_idx(ii),comp_rel_summ(c));
            cpca_struct.Summary(summ_no,c).CompLoadings(ii).idx = opt_comp_load_idx(ii);
        end
        
        if summaryType == 1
            rel_pred_idx = find(cpca_struct.RelPredLoadMask(:,comp_rel_summ(c)));
        elseif summaryType == 2
            rel_pred_idx = find(cpca_struct.RelPredLoadMask_corrected(:,comp_rel_summ(c)));
        else
            error('This should never happen!');
        end
        

        for ii = 1:length(rel_pred_idx)
            cpca_struct.Summary(summ_no,c).PredLoadings(ii).Name = params.Glabel{rel_pred_idx(ii)};
            cpca_struct.Summary(summ_no,c).PredLoadings(ii).Value = cpca_struct.mean_PCorr(rel_pred_idx(ii),comp_rel_summ(c));
            cpca_struct.Summary(summ_no,c).PredLoadings(ii).pValue = cpca_struct.pVales(rel_pred_idx(ii),comp_rel_summ(c));
            cpca_struct.Summary(summ_no,c).PredLoadings(ii).corrected_pValue = cpca_struct.pVales_corrected(rel_pred_idx(ii),comp_rel_summ(c));
            cpca_struct.Summary(summ_no,c).PredLoadings(ii).PredLoadRelProp = cpca_struct.mean_PCorr_reliable(rel_pred_idx(ii),comp_rel_summ(c));
            cpca_struct.Summary(summ_no,c).PredLoadings(ii).idx = rel_pred_idx(ii);
        end

    end
    
    if ~params.default_flag
        Add_more_summary_flag = input('Would you like to select different set of component loadings? 1 - Yes, 0 - No: ');
    else
        Add_more_summary_flag = false;
    end
    summ_no = summ_no + 1;

end
end

function [reportObjects] = reportsAddSummary(reportObjects, inStruct)

%% Add variance table

import mlreportgen.dom.*;
import mlreportgen.report.*;

%% Intialize stuff

Z_label = inStruct.Z_label;
G_label = inStruct.G_label;
Sub_label = inStruct.Sub_label;

fullZ_beh_names = {'fullZ','beh'};
pca_cpca_names = {'PCA','CPCA'};

cutoff_string = {'predrel_cutoff','predrel_cutoff_rel'};
summary_type_string = {'based on p<=0.005','based on multiple comparison corrected p <= 0.05'};

%%

for fullz_beh = 1:2
    for pca_cpca = 1:2

        H = Heading1(Text('Summary: '));
        add(reportObjects.rptObj(fullz_beh,pca_cpca), H);

        for zi = 1:length(Z_label)
            for gi = 1:length(G_label)
                for si = 1:length(Sub_label)

                    result_struct = eval(sprintf("inStruct.%s_%s_struct(zi,gi,si)",pca_cpca_names{pca_cpca},fullZ_beh_names{fullz_beh}));
                    predrel_cutoff = eval(sprintf('result_struct.%s',cutoff_string{fullz_beh}));

                    H = Heading1(Text(sprintf('Z: %s, G: %s, Sub: %s',...
                        Z_label{zi},G_label{gi},Sub_label{si})));
                    add(reportObjects.rptObj(fullz_beh,pca_cpca), H);

                    if ~isempty(fieldnames(result_struct.Summary))
                        for summ_no = 1:size(result_struct.Summary,1)

                            H = Heading2(Text(sprintf('Summary Number : %d, %s',summ_no,summary_type_string{result_struct.Summary(summ_no,1).SummaryType})));
                            add(reportObjects.rptObj(fullz_beh,pca_cpca), H);

                            for comp = 1:size(result_struct.Summary,2)
                                H = Heading3(Text(sprintf('Comp %d (Summary %d)',...
                                    result_struct.Summary(summ_no,comp).CompNo, summ_no)));
                                add(reportObjects.rptObj(fullz_beh,pca_cpca), H);

                                CompLoadings = result_struct.Summary(summ_no,comp).CompLoadings;
                                H = Heading4(Text(sprintf('Dominant Component Loadings (r), %d selected', length(CompLoadings))));
                                add(reportObjects.rptObj(fullz_beh,pca_cpca), H);
                                for jj = 1:length(CompLoadings)
                                    T = Text(sprintf('%s (%.2f) \n', CompLoadings(jj).Name, CompLoadings(jj).Value)); 
                                    add(reportObjects.rptObj(fullz_beh,pca_cpca), T);
                                end
                                
                                if result_struct.Summary(summ_no,1).SummaryType == 1
                                    H = Heading4(Text(sprintf('Reliable Predictor loadings (p<=0.005, PL >= %.4f)', predrel_cutoff)));
                                    add(reportObjects.rptObj(fullz_beh,pca_cpca), H);

                                    PredLoadings = result_struct.Summary(summ_no,comp).PredLoadings;
                                    for jj = 1:length(PredLoadings)
                                        T = Text(sprintf('%s (PL = %.2f, p = %.3f) \n', PredLoadings(jj).Name, ...
                                            PredLoadings(jj).Value, PredLoadings(jj).pValue));
                                        add(reportObjects.rptObj(fullz_beh,pca_cpca), T);
                                    end
                                elseif result_struct.Summary(summ_no,1).SummaryType == 2
                                    H = Heading4(Text(sprintf('Reliable Predictor loadings (corrected p<=0.05)')));
                                    add(reportObjects.rptObj(fullz_beh,pca_cpca), H);

                                    PredLoadings = result_struct.Summary(summ_no,comp).PredLoadings;
                                    for jj = 1:length(PredLoadings)
                                        T = Text(sprintf('%s (PL = %.2f, PLRP = %.2f, p = %.3f) \n', PredLoadings(jj).Name, ...
                                            PredLoadings(jj).Value, PredLoadings(jj).PredLoadRelProp, PredLoadings(jj).corrected_pValue));
                                        add(reportObjects.rptObj(fullz_beh,pca_cpca), T);
                                    end
                                end

                                br = PageBreak();
                                add(reportObjects.rptObj(fullz_beh,pca_cpca), br);
                            end
                        end
                    else
                        T = Text(sprintf('No reliable components!!'));
                        add(reportObjects.rptObj(fullz_beh,pca_cpca), T);
                    end
                end
            end
        end

    end
end

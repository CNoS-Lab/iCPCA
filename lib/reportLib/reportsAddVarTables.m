function [reportObjects] = reportsAddVarTables(reportObjects, inStruct)

%% Add variance table

import mlreportgen.dom.*;
import mlreportgen.report.*;

%% Intialize stuff

Z_label = inStruct.Z_label;
G_label = inStruct.G_label;
Sub_label = inStruct.Sub_label;

fullZ_beh_names = {'fullZ','beh'};
pca_cpca_names = {'PCA','CPCA'};

%%

for fullz_beh = 1:2
    for pca_cpca = 1:2

        result_struct = eval(sprintf("inStruct.%s_%s_struct",pca_cpca_names{pca_cpca},fullZ_beh_names{fullz_beh}));

        for zi = 1:length(Z_label)
            for gi = 1:length(G_label)
                for si = 1:length(Sub_label)

                    paramsTemp = result_struct.params;
                    row_idx = (zi-1)*length(Sub_label)*length(G_label)+(gi-1)*length(Sub_label)+si;

                    table_mat{row_idx,1} = categorical(string([Z_label{zi}, sprintf(' (%d)',paramsTemp.Zdim)]));
                    VarNames{1} = 'Z (nZ)';
                    table_mat{row_idx,2} = categorical(string([G_label{gi}, sprintf(' (%d)',paramsTemp.Gdim)]));
                    VarNames{2} = 'G (nG)';
                    table_mat{row_idx,3} = categorical(string([Sub_label{si}(1), sprintf(' (%d)',paramsTemp.n_sub)]));
                    VarNames{3} = 'Sub (n)';
                    table_mat{row_idx,4} = categorical(string(sprintf('%.2f',paramsTemp.n_sub/paramsTemp.Gdim)));
                    VarNames{4} = 'nSub/nG';

                    varTable = cat(1,result_struct(zi,gi,si).varianceTable);

                    table_mat{row_idx,5} = categorical(string(mean([varTable(:).perVarS],'omitnan')));
                    VarNames{5} = '%Var';

                    compVar = mean(cat(1,varTable(:).perVarScoresCompWiseBetween),1,'omitnan');
                    table_mat{row_idx,6} = categorical(string(compVar(1)));
                    VarNames{6} = 'C1 %Var';

                    if length(compVar) > 1
                        table_mat{row_idx,7} = categorical(string(compVar(2)));
                        VarNames{7} = 'C2 %Var';
                    else
                        table_mat{row_idx,7} = categorical('NA');
                        VarNames{7} = 'C2 %Var';
                    end

                end
            end
        end

        T = cell2table(table_mat);
        T.Properties.VariableNames = VarNames;

        H = Heading4(Text('Summary table for default settings')); 
        add(reportObjects.rptObj(fullz_beh,pca_cpca), H);

        Tbl = MATLABTable(T);
        Tbl.TableEntriesStyle = {Height('20pt')};
        add(reportObjects.rptObj(fullz_beh,pca_cpca), Tbl);

    end
end

end

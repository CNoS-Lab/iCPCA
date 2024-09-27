function reportObjects = reportsInit(savePath, dataSetName, codeVersion, codeDate, r)

import mlreportgen.dom.*;
import mlreportgen.report.*;

fullZ_beh_names = {'fullZ','Iterative'};
pca_cpca_names = {'PCA','CPCA'};

for fullz_beh = 1:2
    for pca_cpca = 1:2
        reportObjects.Name{fullz_beh,pca_cpca} = fullfile(savePath, sprintf('%s_dataset_%s_%s_analysis_%s_r=%0.3f',...
            dataSetName,fullZ_beh_names{fullz_beh},pca_cpca_names{pca_cpca},codeVersion, r));
        reportObjects.Title{fullz_beh,pca_cpca} = sprintf('%s dataset: %s %s analysis %s, r=%.3f',...
            dataSetName,fullZ_beh_names{fullz_beh},pca_cpca_names{pca_cpca},codeVersion, r);

        reportObjects.rptObj(fullz_beh,pca_cpca) = Report(reportObjects.Name{fullz_beh,pca_cpca},'pdf');
        
        tp = TitlePage;
        tp.Subtitle = sprintf('%s: %s',codeVersion, codeDate), ;
        tp.Title = reportObjects.Title{fullz_beh,pca_cpca};
        add(reportObjects.rptObj(fullz_beh,pca_cpca),tp);

        toc = TableOfContents;
        toc.NumberOfLevels = 4;
        add(reportObjects.rptObj(fullz_beh,pca_cpca),toc);

    end
end

end
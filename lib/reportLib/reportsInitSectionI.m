function reportObjects = reportsInitSectionI(reportObjects, idxNos, TitleStr)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2
        reportObjects.SectionsI(idxNos(1), idxNos(2)).rptObj(fullz_beh, pca_cpca) = Section;
        reportObjects.SectionsI(idxNos(1), idxNos(2)).rptObj(fullz_beh, pca_cpca).Title = TitleStr;
    end
end

end
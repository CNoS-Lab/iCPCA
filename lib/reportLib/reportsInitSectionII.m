function reportObjects = reportsInitSectionII(reportObjects, IdxNos, TitleStr)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2
        reportObjects.SectionsII(IdxNos(1), IdxNos(2), IdxNos(3)).rptObj(fullz_beh, pca_cpca) = Section;
        reportObjects.SectionsII(IdxNos(1), IdxNos(2), IdxNos(3)).rptObj(fullz_beh, pca_cpca).Title = TitleStr;
    end
end

end
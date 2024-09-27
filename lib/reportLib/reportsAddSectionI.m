function reportObjects = reportsAddSectionI(reportObjects, IdxNos)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2

        add(reportObjects.Chapters(IdxNos(1)).rptObj(fullz_beh, pca_cpca), ...
            reportObjects.SectionsI(IdxNos(1), IdxNos(2)).rptObj(fullz_beh, pca_cpca));

    end
end

end
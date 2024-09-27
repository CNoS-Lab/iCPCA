function reportObjects = reportsAddSectionII(reportObjects, IdxNos)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2

        add(reportObjects.SectionsI(IdxNos(1), IdxNos(2)).rptObj(fullz_beh, pca_cpca), ...
            reportObjects.SectionsII(IdxNos(1), IdxNos(2), IdxNos(3)).rptObj(fullz_beh, pca_cpca));

        br = PageBreak();
        add(reportObjects.SectionsI(IdxNos(1), IdxNos(2)).rptObj(fullz_beh, pca_cpca), br);

    end
end

end
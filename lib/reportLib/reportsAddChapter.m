function reportObjects = reportsAddChapter(reportObjects, IdxNos)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2

        add(reportObjects.rptObj(fullz_beh,pca_cpca), ...
            reportObjects.Chapters(IdxNos(1)).rptObj(fullz_beh,pca_cpca));

    end
end

end
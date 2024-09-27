function reportObjects = reportsInitChapter(reportObjects, chNo, TitleStr)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2
        reportObjects.Chapters(chNo).rptObj(fullz_beh,pca_cpca) = Chapter;
        reportObjects.Chapters(chNo).rptObj(fullz_beh,pca_cpca).Title = TitleStr;
    end
end

end
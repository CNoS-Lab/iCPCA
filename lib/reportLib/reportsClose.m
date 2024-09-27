function reportObjects = reportsClose(reportObjects)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2

        close(reportObjects.rptObj(fullz_beh,pca_cpca)); 
        rptview(reportObjects.rptObj(fullz_beh,pca_cpca)); 
        
    end
end

end
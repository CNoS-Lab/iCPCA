function rptObj = reportsAddObj(rptObj, Obj)

import mlreportgen.dom.*;
import mlreportgen.report.*;

for fullz_beh = 1:2
    for pca_cpca = 1:2

        add(rptObj(fullz_beh,pca_cpca), Obj);

    end
end

end
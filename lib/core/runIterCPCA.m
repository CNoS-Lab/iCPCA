function [CPCA_struct] = runIterCPCA(Z, G, nComp, params)

nIter = params.nIter;
for ii = 1:nIter
   cv_idx(:,ii) = crossvalind('KFold',size(Z,1),2); 
end

Z = (Z - mean(Z,1,"omitnan"))./std(Z,[],1,"omitnan"); % Z-score the data 

parfor ii = 1:nIter
    for k = 1:2
        if strcmp(params.cpca_pca_flag, 'pca')
            temp = runCPCA(Z(cv_idx(:,ii) == k,:),1,1,2,0);
            temp.G_orig = G(cv_idx(:,ii) == k,:);
            allIterations(ii,k) = temp;
        elseif strcmp(params.cpca_pca_flag, 'cpca')
            temp = runCPCA(Z(cv_idx(:,ii) == k,:),G(cv_idx(:,ii) == k,:),1,2,0);
            temp.G_orig = G(cv_idx(:,ii) == k,:);
            allIterations(ii,k) = temp;
        end
    end
end

S_big = cat(1, allIterations(:).S);

if isempty(nComp)
    var_struct = runCPCA(S_big,1,1,rank(S_big),0);
    plotScree(var_struct);
    nComp = input('Behavioral CPCA (with iterations) - Enter no. of comps: ');
    CPCA_struct = runCPCA(S_big,1,1,nComp,params.varimaxFlag);
else
    CPCA_struct = runCPCA(S_big,1,1,nComp,params.varimaxFlag);
end
CPCA_struct.nComp = nComp;
CPCA_struct.params = params;

CPCA_struct.PCorr = NaN(size(G,2), nComp, nIter, 2);
CPCA_struct.PCorr_pVal = NaN(size(G,2), nComp, nIter, 2);
CPCA_struct.PCorr_reliable = NaN(size(G,2), nComp, nIter);

first_idx = 1;
for k = 1:2
    for ii = 1:nIter
        last_idx = first_idx + size(allIterations(ii,k).Z,1) - 1;
        CPCA_struct.scores_NU_GH_cell{ii,k} = CPCA_struct.scores_NU_GH(first_idx:last_idx,:);
        [CPCA_struct.PCorr(:,:,ii,k), CPCA_struct.PCorr_pVal(:,:,ii,k)] = corr(allIterations(ii,k).G_orig, CPCA_struct.scores_NU_GH_cell{ii,k});
        first_idx = last_idx + 1;
        
        var_struct.Z = allIterations(ii,k).Z;
        var_struct.S = allIterations(ii,k).S;
        var_struct.E = allIterations(ii,k).E;
        var_struct.scores_NU_GH = CPCA_struct.scores_NU_GH_cell{ii,k};
        var_struct.loadings_VD_sqrtN_GH = CPCA_struct.loadings_VD_sqrtN_GH;

        CPCA_struct.varianceTable(ii,k) = getVarianceTable(var_struct);

        if k == 2
            CPCA_struct.PCorr_reliable(:,:,ii) = abs(CPCA_struct.PCorr(:,:,ii,1)) >= params.r & ...
                abs(CPCA_struct.PCorr(:,:,ii,2)) >= params.r;
        end
    end
end

CPCA_struct.mean_PCorr = mean(reshape(CPCA_struct.PCorr, [size(G,2), nComp, nIter*2]),3,"omitnan");
CPCA_struct.mean_PCorr_reliable = mean(CPCA_struct.PCorr_reliable,3,"omitnan");
CPCA_struct.allIterations = allIterations;

CPCA_struct.cv_idx = cv_idx;
CPCA_struct.iterScores = NaN(size(Z,1),nComp,nIter);

for ii = 1:nIter
    CPCA_struct.iterScores(cv_idx(:,ii) == 1,:,ii) = CPCA_struct.scores_NU_GH_cell{ii,1};
    CPCA_struct.iterScores(cv_idx(:,ii) == 2,:,ii) = CPCA_struct.scores_NU_GH_cell{ii,2};
end

CPCA_struct.meanScores = mean(CPCA_struct.iterScores, 3, "omitnan");

end

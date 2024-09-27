function Res = SelectCompLoadingsBig(cpca_struct, params)

Z_dim = size(cpca_struct.Z,2);
S = cpca_struct.S;
G = params.G_orig;

for jj = 1:Z_dim

    Rej_feat = jj;
    Sel_feat = setdiff(1:Z_dim, Rej_feat);

    N = size(S,1);
    G_CPCA = runCpcaStep1(S(:,Sel_feat), S(:,Rej_feat), 1);
    Res(jj).ResidueS = G_CPCA.E;
    Res(jj).ResidueScores = sqrt(N-1) .* Res(jj).ResidueS * cpca_struct.VGH(Sel_feat,:) / cpca_struct.DGH;

    if params.varimaxFlag
        Res(jj).ResidueScoresRot = Res(jj).ResidueScores(:,1:cpca_struct.nComp)*cpca_struct.T;
    else
        Res(jj).ResidueScoresRot = Res(jj).ResidueScores(:,1:cpca_struct.nComp);
    end

    % Break-up the big component score matrix into individual chunks
    % and estimate predictor correlations

    nComp = cpca_struct.nComp;
    nIter = params.nIter;

    allIterations = cpca_struct.allIterations;
    temp_PCorr = NaN(size(G,2), nComp, nIter, 2);
    temp_PCorr_pVal = NaN(size(G,2), nComp, nIter, 2);
    temp_PCorr_reliable = NaN(size(G,2), nComp, nIter);

    first_idx = 1;
    for k = 1:2
        for ii = 1:nIter
            last_idx = first_idx + size(allIterations(ii,k).Z,1) - 1;
            temp_scores = Res(jj).ResidueScoresRot(first_idx:last_idx,:);
            [temp_PCorr(:,:,ii,k), temp_PCorr_pVal(:,:,ii,k)] = corr(allIterations(ii,k).G_orig, temp_scores);
            first_idx = last_idx + 1;

            if k == 2
                temp_PCorr_reliable(:,:,ii) = temp_PCorr_pVal(:,:,ii,1) <= 0.05 & ...
                    temp_PCorr_pVal(:,:,ii,2) <= 0.05;
            end
        end
    end

    Res(jj).mean_PCorr_reliable = mean(temp_PCorr_reliable,3,"omitnan");
    Res(jj).ResiduePCorr = mean(reshape(temp_PCorr, [size(G,2), nComp, nIter*2]),3,"omitnan");
    
    Res(jj).meanPCorr_rel = sum(abs(Res(jj).mean_PCorr_reliable).*cpca_struct.RelPredLoadMask, 1)./sum(cpca_struct.RelPredLoadMask, 1);
    Res(jj).meanPCorr_rel_corrected = sum(abs(Res(jj).mean_PCorr_reliable).*cpca_struct.RelPredLoadMask_corrected, 1)./sum(cpca_struct.RelPredLoadMask_corrected, 1);

end
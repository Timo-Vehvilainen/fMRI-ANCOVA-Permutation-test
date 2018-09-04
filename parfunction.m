function [pval tstat_unpermuted] = parfunction(memMaps_c1, memMaps_c2, x, y, z)
%PARFUNCTION 
    % Used for performing the statistical analysis (Both the ANCOVA and the
    % permutation testing) inside a parallizable for-loop

    %Number of subjects per group
    Nsubs = 29; %Number of total subjects
    G1_Nsubs=15; % Number of subjects for group 1
    G2_Nsubs=14; % Number of subjects for group 2
    T=200; % Number of time points

    %lists of indices for both groups
    G1_ids = 1:G1_Nsubs;
    G2_ids = (G1_Nsubs + 1):Nsubs;  

    %Constructing indicator-matrices for each group, 
    %with the appropriate indices flagged
    G1 = zeros(Nsubs);
    G1(G1_ids, G1_ids) = 1;
    G2 = zeros(Nsubs);
    G2(G2_ids, G2_ids) = 1;
    
    %indices for upper triangular matrix (used for extracting the data)
    UT_ids = find(triu(ones(Nsubs), 1));
    
    %First for C1
    g12c1_v = ...
        squeeze(memMaps_c1.memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc(x,y,z,:));
    g12c1 = zeros(Nsubs);
    %change from vector form to upper triangular matrix
    g12c1(UT_ids) = g12c1_v;  
    %mirror the lower triangular aswell, and insert the diagonal
    g12c1 = g12c1 + g12c1' + eye(size(g12c1)); 

    %Then for C2
    g12c2_v = ...
        squeeze(memMaps_c2.memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc(x,y,z,:));
    g12c2 = zeros(Nsubs);
    g12c2(UT_ids) = g12c2_v;
    g12c2 = g12c2+g12c2'+eye(size(g12c2));

    %Starting the ANCOVA step
    %Extract the lower triangular of the original matrices in vector form
    LT_ids = find(tril(ones(Nsubs), -1));
    C1 = atanh(g12c1(LT_ids));
    C2 = atanh(g12c2(LT_ids));

    %Regress C1 out of C2 (the ANCOVA step). 
    %The leftover residual variation is in R.
    [Beta, ~, R] = regress(C2,[C1 ones(size(C1))]);

    %Reshape R into a symmetric matrix in preparation for permutations
    Data = ones(Nsubs);
    Data(LT_ids) = R;
    Data = Data + Data';    %Adding also the top triangle before permuting
    %End of ANCOVA step

    %Starting the PERMUTATION TEST
    %Extract groups 1 and 2 from R
    g1_LT_ids = find(tril(G1, -1));
    g1 = Data(g1_LT_ids);
    g2_LT_ids = find(tril(G2, -1));
    g2 = Data(g2_LT_ids);

    %Record the t-stat value for the unpermuted voxel
    [~, ~, ~, STATS] = ttest2(g1, g2);
    tstat_unpermuted = STATS.tstat;

    Niter = 5000;

    %initialize matrix to hold the permutations and the permuted t-stats
    perms = zeros(5000, size(Data, 1));
    tstats = zeros(Niter, 1);

        %Generate random permutations
        for N = 1:Niter
            perms(N, :) = randperm(size(Data, 1))';
        end

        for N = 1:Niter
            %Implement permutation on R
            fake_Data = Data(perms(N, :), perms(N, :));

            %Extract groups 1 and 2
            fake_g1 = fake_Data(g1_LT_ids);
            fake_g2 = fake_Data(g2_LT_ids);

            %Calculate and record t-stat for each permutation
            [~, ~, ~, STATS] = ttest2(fake_g2, fake_g1);  
            tstats(N, :) = STATS.tstat;

        end

    %calculate p-value for unpermuted t-stat by interpolating the cdf of tstats
    [F XI] = ksdensity(tstats, 'function', 'cdf', 'npoints', 200);
    pval_left = interp1([-1 XI 1], [0 F 1], tstat_unpermuted);
    pval = 1-pval_left;

    %End of PERMUTATION TEST

end


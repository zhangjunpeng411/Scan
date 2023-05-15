%%%%% Running GenMiR++ in K562 dataset %%%%%
%% No prior knowledge
t1 = clock;
core_number = 32;
p_value_cutoff = 0.05;
nsamples = size(Z, 2);
parpool(core_number);
parfor i = 1:size(Z, 2)
    Z_update = Z;
    X_update = X;
    Z_update(:, i) = [];
    X_update(:, i) = [];
    Score_single_null{i} = GenMiR_generic(Z_update, X_update, C_null);
end
Score_all_null = GenMiR_generic(Z, X, C_null);
for i = 1:size(Z, 2)
    res_single_null{i} = nsamples * Score_all_null - (nsamples - 1) * Score_single_null{i}; 
    res_single_null{i} = 1- normcdf(abs(matrixzcore(res_single_null{i})));
    res_single_null_pvalue{i} = res_single_null{i} < p_value_cutoff;
end

delete(gcp('nocreate'))
t2 = clock;
time_null = etime(t2,t1);

%% Extract results
filepattern1 = 'Scan.interp_K562_NULL%-d.csv'; 
for iter = 1:numel(res_single_null_pvalue)
    csvwrite(sprintf(filepattern1, iter), res_single_null_pvalue{1, iter});
end

%%%%% Running GenMiR++ in BRCA dataset %%%%%
%% No prior knowledge
t1 = clock;
core_number = 32;
p_value_cutoff = 0.05;
nsamples = size(Z, 2);
parpool(core_number);
parfor i = 1:size(Z, 2)
    Z_update = Z;
    X_update = X;
    Z_update(:, i) = [];
    X_update(:, i) = [];
    Score_single_null{i} = GenMiR_generic(Z_update, X_update, C_null);
end
Score_all_null = GenMiR_generic(Z, X, C_null);
for i = 1:size(Z, 2)
    res_single_null{i} = nsamples * Score_all_null - (nsamples - 1) * Score_single_null{i}; 
    res_single_null{i} = 1- normcdf(abs(matrixzcore(res_single_null{i})));
    res_single_null_pvalue{i} = res_single_null{i} < p_value_cutoff;
end

delete(gcp('nocreate'))
t2 = clock;
time_null = etime(t2,t1);

%% Extract results
filepattern2 = 'Scan.interp_BRCA_NULL%-d.csv'; 
for iter = 1:numel(res_single_null_pvalue)
    csvwrite(sprintf(filepattern1, iter), res_single_null_pvalue{1, iter});
end


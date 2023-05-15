%%%%% Running GenMiR++ in K562 dataset %%%%%
%% No prior knowledge
t1 = clock;
core_number = 32;
parpool(core_number);
parfor i = 1:size(Z, 2)
    Z_update = Z;
    X_update = X;
    Z_update(:, i) = [];
    X_update(:, i) = [];
    Score_single_null{i} = GenMiR_SSN(Z_update, X_update, C_null, 0.05);
end
Score_all_null = GenMiR_SSN(Z, X, C_null, 0.05);
for i = 1:size(Z, 2)
    res_single_null{i} = abs(Score_single_null{i} - Score_all_null);
end

delete(gcp('nocreate'))
t2 = clock;
time_null = etime(t2,t1);

%% Extract results
filepattern1 = 'Scan.perturb_K562_NULL%-d.csv'; 
for iter = 1:numel(res_single_null)
    csvwrite(sprintf(filepattern1, iter), res_single_null{1, iter});
end

%%%%% Running GenMiR++ in BRCA dataset %%%%%
%% No prior knowledge
t1 = clock;
core_number = 32;
parpool(core_number);
parfor i = 1:size(Z, 2)
    Z_update = Z;
    X_update = X;
    Z_update(:, i) = [];
    X_update(:, i) = [];
    Score_single_null{i} = GenMiR_SSN(Z_update, X_update, C_null, 0.05);
end
Score_all_null = GenMiR_SSN(Z, X, C_null, 0.05);
for i = 1:size(Z, 2)
    res_single_null{i} = abs(Score_single_null{i} - Score_all_null);
end

delete(gcp('nocreate'))
t2 = clock;
time_null = etime(t2,t1);

%% Extract results
filepattern2 = 'Scan.perturb_BRCA_NULL%-d.csv'; 
for iter = 1:numel(res_single_null)
    csvwrite(sprintf(filepattern1, iter), res_single_null{1, iter});
end

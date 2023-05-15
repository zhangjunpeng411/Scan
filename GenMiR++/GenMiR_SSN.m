function res = GenMiR_SSN(Z, X, C, p_value_cutoff)

res_origin = GenMiR_generic(Z, X, C);
C(C == 0) = NaN; 
res_origin = res_origin .* C;
res_zscore = 1 - normcdf(matrixzcore(res_origin));
res = res_zscore < p_value_cutoff;

end

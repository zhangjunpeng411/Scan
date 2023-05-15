function res = GenMiR_generic(Z, X, C)

res = zeros(size(C));
[row col] = find(C);
Score = GenMiR(Z, X, C);
for i=1:length(row)
    res(row(i),col(i)) = Score(i); %The matrix with rows mRNAs, and columns miRNAs    
end

end

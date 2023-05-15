function  res = matrixzcore(mat)
    mean_value = nanmean(mat(:));
    std_value = nanstd(mat(:));
    res = (mat -mean_value)/std_value;
end
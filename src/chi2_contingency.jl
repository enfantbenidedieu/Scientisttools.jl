using Distributions

function chi2_contingency(X)
    observed = Matrix{Int}(X);
    l, c, total = size(observed)[1],size(observed)[2],sum(observed)
    rowsum, colsum = sum(observed,dims=2),sum(observed,dims=1)
    expected_freq = (reshape(rowsum,(l,1)) * reshape(colsum,(1,c)))/total
    statistic, dof = sum((observed - expected_freq).^2 ./ expected_freq), (l-1)*(c-1)
    p_value = 1 - cdf(Chisq(dof),statistic)
    return (; :statistic => statistic, :pvalue => p_value, :dof => dof, :expected_freq => expected_freq)
end;
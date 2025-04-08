# Coefficient V de Cramer
function cramers_v(phi2,l,c)
    return sqrt(phi2/min(l-1,c-1))
end;

#coefficient T de Tschuprow
function tschuprow_t(phi2,l,c)
    return sqrt(phi2/sqrt((l-1)*(c-1)))
end;

#coefficient C de Pearson
function pearson_c(phi2)
    return sqrt(phi2/(phi2+1))
end;

"""
    association_measure(X)

Compute association measure on continency table

...
# Description



...
"""
function association_measure(X)
    X = Matrix{Float64}(X);
    l,c,n = size(X)[1],size(X)[2],sum(X)
    chi2,p_value,dof,expected_freq = chi2_contingency(X)
    phi2 = chi2/n
    coef_v,coef_t,coef_c = cramers_v(phi2,l,c),tschuprow_t(phi2,l,c),pearson_c(phi2)
    #
    res = hcat(DataFrame(Mesure=["Coefficient"]),
               DataFrame([chi2 p_value phi2 coef_v coef_t coef_c],
                        ["chi2","p-value","phi2","V Cramer","T Tschuprow","C Pearson"]))
    return res
end;


function g_test(X)
    observed = Matrix{Int}(X);
    l,c,total = size(observed)[1],size(observed)[2],sum(observed)
    rowsum, colsum = sum(observed,dims=2),sum(observed,dims=1)
    expected_freq = (reshape(rowsum,(l,1)) * reshape(colsum,(1,c)))/total
    iar = observed ./ expected_freq
    statistic, dof = 2 * sum(observed .* log.(iar)),  (l-1)*(c-1);
    p_value = 1 - cdf(Chisq(dof),statistic)
    return  (; :statistic => statistic, :pvalue => p_value, :dof => dof)
end;

using XLSX,DataFrames
data,columns = XLSX.readtable("./data/children.xlsx","Feuil1");
donnee = DataFrame(data,columns);
association = association_measure(donnee[!,2:end]);
chi2_test = chi2_contingency(donnee[!,2:end]);
gtest = g_test(donnee[!,2:end]);
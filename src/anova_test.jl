"""
    anova_test(categories,values)

Analysis of variance

...
# Description
Compute analysis of variance between categories variables and quantitative variables.

# Arguments
- `categories`: a vector of categories data
- `values`:  a vector of quantitative data

# Returns
A DataFrame containing inter variance, intra variance, square correlation ratio, F-stats and p-value associated.

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function anova_test(categories::AbstractVector,
                    values::AbstractVector)::AbstractDataDFrame
    u = sort(unique(categories))
    q = length(u)
    var_intra, var_inter = 0.0,0.0
    for cat in u
        subgroup=values[findall(x->x==cat,vec(categories))]
        var_intra=var_intra+sum((subgroup .- mean(subgroup)) .^2)
        var_inter=var_inter+length(subgroup)*((mean(subgroup)-mean(values)).^2)
    end;
    eta2 = var_inter/(var_inter + var_intra)
    ddl1, ddl2 = q - 1, length(values) - q
    f_stats = (ddl2/ddl1)*(var_inter/var_intra)
    f_qt = quantile(FDist(ddl1,ddl2),0.95) #quantile
    p_value = 1 - cdf(FDist(ddl1,ddl2),f_stats) #pvalue

    # Store all results and convert to R2
    stats = [var_inter var_intra eta2 f_stats f_qt p_value]
    return DataFrame(stats,["Inter Var","Intra Var.","R2","F-statistic","F-value","p-value"])
end;
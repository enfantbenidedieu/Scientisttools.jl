"""
    contdesc(X,y;...)

Continuous variables description

...
# Description
Description continuous by quantitative variables

# Required arguments
- `X`: a DataFrame of shape (`n_samples`, `n_columns`) of quantitative
- `y`: a vector of factor coordinate of length `n_samples`.

# Optional (default) arguments
- `w`: a optional individuals weights (default = nothing)
- `proba`: the significance threshold considered to characterized the category (by default 0.05)

# Returns
A DataFrame (n_columns, 2) containing weighted correlation and the pvalue associated :
- `statistic`: weighted Pearson product-moment correlation coefficient 
- `p-value`: the p-value associated

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function contdesc(X::AbstractDataFrame,
                  y::AbstractVector;
                  w=nothing,
                  proba::Float64=0.05)::AbstractDataFrame

    if w === nothing
        w = ones(nrow(X))/nrow(X)
    else
        w = w / sum(w)
    end;

    # Fill NA with mean
    X = recodecont(X).Xcod

    function cor_test(x,y,w,col)
        return hcat(DataFrame("Variables" => col),wcortest(x,y,w=w)[!,["statistic","p-value"]])
    end;  

    value = reduce(vcat,[cor_test(vec(X[!,col]),vec(y),w,col) for col in names(X)])
    # Filter
    value = filter(row -> row["p-value"] < proba, value)
    # Sort by statistic
    value = sort(value,Symbol("statistic"), rev=true)
    # Rename
    rename!(value, :statistic => :correlation) 
    return value
end;
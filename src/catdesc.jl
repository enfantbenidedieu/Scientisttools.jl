"""
    catdesc(X,y;...)

Categories description

...
# Description
Description of the categories of one factor by categorical variables;

# Required arguments
- `X`: a DataFrame of shape (`n_samples`, `n_columns`) of qualitative variables
- `y`: a vector of factor coordinate of length `n_samples`.

# Optional (default) arguments
- `proba`: the significance threshold considered to characterized the category (by default 0.05)

# Retuns
a NamedTuple containing:
- `quanti`: a DataFrame of square correlation ratio and p-value associated
- `category`: a DataFrame of estimated coefficients and p-value associated.

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function catdesc(X::AbstractDataFrame,
                 y::AbstractVector;
                 proba::Float64=0.05)::NamedTuple

    # Test if X is DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("'X' must be a DataFrame.")) end;

    # Test if y is vector
    if !isa(y,AbstractVector) throw(ArgumentError("'y' must be a vector.")) end;

    # Analysis of variables
    function model(y,x,col)
        return hcat(DataFrame("variables" => col),anova_test(x,y)[!,["R2","p-value"]])
    end;
    quali = reduce(vcat,[model(vec(y),vec(X[!,col]),col) for col in names(X)])
    # Filter
    quali = filter(row -> row["p-value"] < proba, quali)
    # Sort by R2
    quali = sort(quali,Symbol("R2"), rev=true)

    function estimate(x,y,col)
        data = DataFrame(X=x,Y=y)
        ols = lm(@formula(Y ~ X), data)
        t_stats = coef(ols)[2]/stderror(ols)[2]
        p_value =  2*(1 - cdf(TDist(dof_residual(ols)),t_stats))
        return hcat(DataFrame("Modalites" => col),DataFrame([coef(ols)[2] p_value],["Estimate","p-value"]))
    end;
    # Categories
    dummies = get_dummies(X,prefix=true,prefix_sep="=")
    category =  reduce(vcat,[estimate(vec(y),vec(dummies[!,col]),col) for col in names(dummies)])

    # Filter
    category = filter(row -> row["p-value"] < proba, category)
    return (; :quali => quali, :category => category)
end;
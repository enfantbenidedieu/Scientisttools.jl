"""
    wcortest(x,y;...)

Weighted Pearson correlation coefficient and p-value for testing non-correlation.

...
# Description
Test for weighted pearson correlation coefficient.

# Required arguments
- `x`: a vector with quantitative values
- `y`: a vector with quantitative values

# Optional (default) arguments
- `w`: a optional individuals weights (default = nothing)

# Returns
a DataFrame containing estimates of the weighted correlation, the degree of freedom and the pvalue associated :
- `statistic`: weighted Pearson product-moment correlation coefficient
- `dof`: degre of freedom   
- `p-value`: the p-value associated

# Examples
```julia-repl
julia> using Scientisttools
julia> x = reduce(vcat,collect.(range(1,10)));
julia> y = [1,2,3,8,7,6,5,8,9,10];
julia> wt = [0,0,0,1,1,1,1,1,0,0]
julia> wcor_test = wcortest(x,y,w=wt);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function wcortest(x::AbstractVector,y::AbstractVector;w=nothing)

    # Test if x is a vector
    if !isa(x,AbstractVector) throw(ArgumentError("'x' must be a vector.")) end;

    # Test if y is a vector
    if !isa(y,AbstractVector) throw(ArgumentError("'y' must be a vector.")) end;

    # Set number of rows
    n_rows = length(x)
    # If w is nothing
    if w === nothing w = ones(n_rows)/n_rows end;
    
    # weighted corr test
    stats = cor(hcat(x,y),weights(w))[1,2]
    t_stats, dof = stats * sqrt((n_rows - 2)/(1 - (stats .^2))), n_rows - 2
    p_value =  2*(1 - cdf(TDist(dof),abs(t_stats)))
    return DataFrame([stats dof p_value],["statistic","dof","p-value"])
end;
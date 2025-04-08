"""
    get_dummies(X;...)

Convert categorical variable into dummy/indicator variables.

...
# Description
Each variable is converted in as many 0/1 variables as there are different values. Columns in the output are each named after a value; if the input is a DataFrame, the name of the original variable is prepended to the value.

# Required arguments
- `X`: a DataFrame of shape (`n_samples`, `n_columns`)

# Optional (default) arguments
- `prefix`: a boolean to append DataFrame columns names (default = true)
- `prefix_sep`: a  string (default "_"). If appending prefix, separator/delimiter to use.

# Examples
```julia-repl
julia> using Scientisttools, DataFrames
julia> s = DataFrame(s=["a","b","c","d"]);
julia> dummies = get_dummies(s);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_dummies(X::AbstractDataFrame;
                    prefix::Bool=false,
                    prefix_sep="_")

    # Test if X is a DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Create a function for a signe vector
    function dummies(x::AbstractVector;lab::AbstractString,prefix::Bool,prefix_sep="_")
        u, m = sort(unique(x)), transpose(indicatormat(x))
        # Update columns names using prefix
        if prefix === true
            v = [lab * prefix_sep for x in 1:length(u)]
            u = map(string,v,u)
        end;
        return DataFrame(Matrix{Int64}(m),u)
    end;

    # Columns names
    res = [dummies(X[!,q],lab=val,prefix=prefix,prefix_sep=prefix_sep) for (q,val) in enumerate(names(X))]
    return reduce(hcat,res)
end;
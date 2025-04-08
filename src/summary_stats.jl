# Load packages
using StatsBase
using Statistics
using DataFrames

function summary_stats(X::AbstractDataFrame)

    # Check if X is a DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Set colnames
    var_names = DataFrame(Variables=[string(k) for k in names(X)])

    # Convert to matrix
    X = Matrix{Float64}(X)

    # Number of columns
    n_vars = nrow(var_names)

    # Statistics
    function statistics(x)
        stats=[minimum(x) mean(x) std(x,corrected=true) std(x,corrected=false) quantile(x,0.25) quantile(x,0.25) quantile(x,0.75) maximum(x)]
        return DataFrame(stats,["min","mean","bias std","std","25%","50%","75%","max"])
    end;
    # Convert to dataframe
    stats = hcat(var_names,reduce(vcat,[statistics(X[:,k]) for k in 1:n_vars]))
    return stats
end;
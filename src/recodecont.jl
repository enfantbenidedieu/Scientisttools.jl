# Load packages
using DataFrames

# Load intern functions
include("missing_imputation.jl")

function recodecont(X::AbstractDataFrame)

    # Check if X is a DataFrame
    if !(X isa AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Remove categorical columns
    quali_col = [names(X)[i] for i in 1:ncol(X) if isa(X[1,i],String)]
    if length(quali_col)>0 X = select(X,Not(quali_col)) end;

    # Test if dataframe is empty
    if ncol(X) === 0 throw(DimensionMismatch("All variables in X must be numeric.")) end;

    # Fill missing with mean
    X = missing_imputation(X)
    
    # Store new dataframe
    Xcod = X
    if ncol(X) == 1
        Z = nothing
        means = X
        stds = nothing 
    else
        means = mean(Matrix{Float64}(Xcod),dims=1)
        stds = std(Matrix{Float64}(Xcod),dims=1,corrected=false)
        Z = (Matrix{Float64}(Xcod) .- means) ./ stds
    end;
    return (; :Z => Z, :means => means, :std => stds, :Xcod => Xcod)
end;
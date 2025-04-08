# Load packages
using DataFrames

# Load intern functions
include("recode_cat_variable.jl")
include("get_dummies.jl")

function recodecat(X::AbstractDataFrame)

    # Check if X is a DataFrame
    if !(X isa AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Remove numeric columns
    quanti_col = [names(X)[i] for i in 1:ncol(X) if isa(X[1,i],Number)]
    if length(quanti_col)>0 X = select(X,Not(quanti_col)) end;

    # Test if dataframe is empty
    if ncol(X) === 0 throw(DimensionMismatch("All variables in X must be either object or category.")) end;

    # revaluate
    X = recode_cat_variable(X)
    dummies = get_dummies(X)
    return (; :X => X, :dummies => dummies)
end;
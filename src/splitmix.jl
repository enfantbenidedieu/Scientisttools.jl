# Load package
using DataFrames

function splitmix(X::AbstractDataFrame)

    # Check if X is a DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;
    
    # Initialisation
    X_quanti, X_quali = nothing, nothing
    
    # Select quantitative columns
    quanti_col = [names(X)[i] for i in 1:ncol(X) if isa(X[1,i],Number)]

    # Select qualitative columns
    quali_col = [names(X)[i] for i in 1:ncol(X) if isa(X[1,i],String)]

    if length(quanti_col)>0 X_quanti = DataFrame(Matrix{Float64}(select(X,quanti_col)),quanti_col) end;
    if length(quali_col)>0 X_quali = DataFrame(Matrix{String}(select(X,quali_col)),quali_col) end;

    return (; :quanti => X_quanti, :quali => X_quali)
end;
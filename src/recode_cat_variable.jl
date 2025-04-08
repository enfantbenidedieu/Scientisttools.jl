# Load packaages
using DataFrames

include("intersect.jl")

function recode_cat_variable(X::AbstractDataFrame)

    # Check if X is a DataFrame
    if !(X isa AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Create a copy of X
    Y = copy(X)
    if ncol(X)>1
        for k in 1:(ncol(X)-1), l in (k+1):ncol(X)
            if any(x -> isa(x,String), X[!,k]) && any(x -> isa(x,String), X[!,l])
                if length(intersect(unique(X[!,k]),unique(X[!,l]))) > 0
                    valuek = Dict(string(i) => names(X)[k] * "_"* string(i) for i in sort(unique(X[!,k])))
                    valuel = Dict(string(i) => names(X)[l] * "_"* string(i) for i in sort(unique(X[!,l])))
                    Y[!,k], Y[!,l] = get.((valuek, ),X[!,k], missing),get.((valuel,),X[!,l], missing)
                end;
            end;
        end;
    end;
    return Y  
end;
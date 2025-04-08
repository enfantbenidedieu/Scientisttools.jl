# Load packages
using DataFrames

function convert_to_float(X::AbstractDataFrame)
    for col in names(X)
        if any(x -> isa(x,String), X[!,col])
            X[!,col] = parse.(Float64,X[!,col])
        elseif (any(x -> isa(x,Any), X[!,col]))
            X[!,col] = convert(Array{Float64,1}, X[!,col])
        end;
    end;
    return X
end; 
using DataFrames
using CategoricalArrays

function convert_to_categorical(X::AbstractDataFrame)
    function to_categorical(x::AbstractVector,lab::AbstractString)
        return DataFrame(Dict(lab => CategoricalArray(x,ordered=true)))
    end;
    return reduce(hcat,[to_categorical(X[!,col],col) for col in names(X)])
end;
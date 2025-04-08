# Load packages
using DataFrames

function freq_table(X::AbstractDataFrame)
    # Check if X is a DataFrame
    if !(X isa AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    function freq(col)
        stats = combine(groupby(X,col),nrow,proprow)
        noms_mod = DataFrame(Modalities=[string(x) for x in stats[:,1]])
        stats = hcat(noms_mod,DataFrame(Matrix(stats[!,2:end]),["Effectifs","Proportions"]))
        sort!(stats,:Modalities)
        return stats
    end;
    return reduce(vcat,[freq(col) for col in names(X)])
end;
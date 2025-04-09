# Load packages
using DataFrames

# Load intern functions
include("splitmix.jl")
include("recodecont.jl")
include("recodecat.jl")

function recodevarfamd(X::AbstractDataFrame)

    # Check if X is a DataFrame
    if !(X isa AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Split X
    X_quanti, X_quali = splitmix(X).quanti, splitmix(X).quali

    # Initialize
    nb_moda, dummies = nothing, nothing

    if X_quanti !== nothing
        X_quanti = recodecont(X_quanti).Xcod
        n1, k1 = nrow(X_quanti), ncol(X_quanti)
    end;

    if X_quali !== nothing
        n2, k2 = nrow(X_quali), ncol(X_quali)
        rec2 = recodecat(X_quali)
        X_quali, dummies = rec2.X, rec2.dummies
        nb_moda = [length(unique(X_quali[!,col])) for col in names(X_quali)]
    end;

    # Collapse result
    if X_quanti !== nothing && X_quali !== nothing X, n, k = hcat(X_quanti,X_quali), n1, k1 + k2 end;

    if X_quanti !== nothing && X_quali === nothing X, n, k, k2 = X_quanti, n1, k1, 0 end;

    if X_quanti === nothing && X_quali !== nothing X, n, k, k1 = X_quali, n2, k2, 0 end;

    # Store all informations
    return (; :X => X, :n => n, :k => k, :k1 => k1, :k2 => k2, :nb_moda => nb_moda, :quanti => X_quanti, :quali => X_quali, :dummies => dummies);
end;

#=
using XLSX
data,columns = XLSX.readtable("./data/decathlon.xlsx","Feuil1");
data2,columns2 = XLSX.readtable("./data/races_canines.xlsx","Feuil1");
donnee = DataFrame(data,columns);
donnee2 = DataFrame(data2,columns2);
rec1 = recodevarfamd(donnee[!,2:end]);
rec2 = recodevarfamd(donnee2[!,2:end]);
=#


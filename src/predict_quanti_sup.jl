# Load packages
using DataFrames

# Load inter functions
include("convert_to_float.jl")
include("summary_stats.jl")
include("descrstatsw.jl")

function predict_quanti_sup(X::AbstractDataFrame,row_weights::AbstractVector,U::AbstractMatrix)

    # Create labels
    col_names = DataFrame(Variables=string.(names(X)))

    # set factor analysis columns
    dim_index = ["Dim." * string(i) for i in 1:size(U)[2]]

    # Convert to float
    X = convert_to_float(X)

    # Summary statistics
    summary = summary_stats(X)
    
    # Convert to matrix
    X = Matrix{Float64}(X)

    # Compute weighted average and weighted standard deviation
    wmeans, wstd = descrstatsw(X,weights(row_weights))

    # Standardisation
    Z = (X .- wmeans) ./ wstd

    # Factor coordinates
    coord = transpose((Z .* row_weights)) * U
    
    # Square distance to origin
    cor_vcs = (Z .^2) .* row_weights
    sqdisto = transpose(cor_vcs) * ones(size(X)[1]) 
    
    # Cos2
    cos2 = (coord .^2) ./ sqdisto

    # Store all informations
    quanti_sup_ = (; :coord => hcat(col_names,DataFrame(coord,dim_index)),
                     :cor => hcat(col_names,DataFrame(coord,dim_index)),
                     :cos2 => hcat(col_names,DataFrame(cos2,dim_index)));

    return  (; :quanti_sup => quanti_sup_, :summary => summary)
end;
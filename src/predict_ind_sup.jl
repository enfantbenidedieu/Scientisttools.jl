# Load packages
using DataFrames

function predict_ind_sup(Z::AbstractMatrix,row_names::AbstractDataFrame,col_weights::AbstractVector,V::AbstractMatrix)

    # set columns
    dim_index = ["Dim." * string(i) for i in 1:size(V)[2]]

    # Factor coordinates
    coord = (Z .* transpose(col_weights)) * V

    # Square distance to origin
    sqdisto = sum((Z .^2) .* transpose(col_weights),dims=2)
    
    # Cos2
    cos2 = (coord .^2) ./sqdisto
    
    # Store all informations
    res = (; :coord => hcat(row_names,DataFrame(coord,dim_index)),
             :cos2 => hcat(row_names,DataFrame(cos2,dim_index)),
             :dist2 => hcat(row_names,DataFrame("Sq. Dist" => vec(sqdisto))))
    return res
end;
using LinearAlgebra

function svd_triplet(X;row_weights=nothing,col_weights=nothing,n_components=nothing)

    # extract dimension
    n_rows, n_cols = size(X)

    # set row weights
    if row_weights === nothing
        row_weights = ones(n_rows)/n_rows
    end;

    # set columns weights
    if col_weights === nothing
        col_weights = ones(n_cols)
    end;

    row_weights,col_weights = vec(row_weights),vec(col_weights)

    if n_components === nothing
        n_components = min(n_rows - 1, n_cols)
    else
        n_components = min(n_components, n_rows - 1, n_cols)
    end;

    # update X
    X = (X .* transpose(sqrt.(col_weights))) .* sqrt.(row_weights)

    # Singular Value Decomposition
    if n_cols < n_rows
        # SVD to X
        U,delta,V = LinearAlgebra.svd(X,full=false)

        # Update 
        if n_components > 1
            # Find sign
            mult = sign.(sum(V,dims=1))

            # Replace signe 0 by 1
            mult[mult .== 0] .= 0

            # Update U and V with mult
            U,V = U .* mult, V .* mult
        end;

        # Update U and V with weights
        U, V = U ./ sqrt.(row_weights), V ./ sqrt.(col_weights)
    else
        # SVD to transpose X
        V, delta, U = LinearAlgebra.svd(transpose(X),full=false)

        # Update 
        if n_components > 1
            # Find sign
            mult = sign.(sum(V,dims=1))

            # Replace signe 0 by 1
            mult[mult .== 0] .= 0

            # Update U and V with mult
            U,V = U .* mult, V .* mult
        end;

        # Update U and V with weights
        U, V = U ./ sqrt.(row_weights), V ./ sqrt.(col_weights)
    end;

    # set number of columns using n_components
    U, V = U[:,1:n_components], V[:,1:n_components]

    # set delta length
    vs = delta[1:min(n_cols, n_rows - 1)]

    # filter to remove value less than 1e-15
    vs_filter = ifelse.(vs[1:n_components] .< 1e-15,1,0)
    
    # select index which respect the criteria
    num = [i for (i, val) in enumerate(vs[1:n_components]) if vs_filter[i]==1]
    vs_num = [val for (i, val) in enumerate(vs[1:n_components]) if vs_filter[i]==1]

    if length(num)==1
        U[:,num], V[:,num] = U[:,num] * vs_num, V[:,num] * vs_num
    elseif  length(num)>1
        U[:,num], V[:,num] = U[:,num] * vs_num, V[:,num] * vs_num
    end;
    
    return (; :vs => vs, :U => U, :V => V)
end;
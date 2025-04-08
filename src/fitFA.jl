"""
    fitFA(Z,max_components;...)

Fit Factor Analysis Method (PCA, CA, MCA)

...
# Required arguments
- `Z`: Matrix of shape (`n_samples`, `n_columns`)
    Training data, where `n_samples` is the number of samples and `n_columns` is the number of columns.
- `max_components`: an integer specified the maximum of dimensions.

# Optional (default) arguments
- `n_components`: an integer specified the number of dimensions kept in the results (by default 5).
- `row_weights`: an optional rows weights (by default, a vector of 1/(number of active rows) for uniform individuals weights), the weights are given only for active rows.
- `row_names`: an optional rows names (by default, a dataframe).
- `col_weights`: an optional columns weights (by default, a vector of 1 for uniform variables weights), the weights are given only for active columns.
- `col_names`: an optional columns names (by default, a dataframe).

# Returns
A NamedTuple with following keys:
- `svd`: namedtuple with results for generalied singular value decomposition (GSVD)
- `eig`: dataFrame containing all the eigenvalues, the difference between each eigenvalue, the percentage of variance and the cumulative percentage of variance.
- `row`: namedtuple of DataFrame containing all the results for the active rows (factor coordinates, square cosinus, relative contributions, additionals informations).
- `col`: namedtuple of DataFrame containing all the results for the active columns (factor coordinates, square cosinus, relative contributions, additionals informations).

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fitFA(Z::AbstractMatrix,
               max_components::Int;
               n_components::Int=5,
               row_weights=nothing,
               row_names=nothing,
               col_weights=nothing,
               col_names=nothing)::NamedTuple

    # Extract size of Z
    n_rows, n_cols = size(Z)

    # Set row and columns weights
    if row_weights === nothing row_weights = ones(n_rows)/n_rows end;
    if col_weights === nothing col_weights = ones(n_cols) end;

    # Set rows and columns names
    if row_names === nothing row_names = DataFrame(Rows=["ind" * string(i) for i in 1:n_rows]) end;
    if col_names === nothing col_names = DataFrame(Columns=["ind" * string(k) for k in 1:n_cols]) end;
    
    # set columns for factor coordintes
    dim_index = ["Dim." * string(i) for i in 1:n_components]

    #########################################################################################################
    ## Rows informations : weights, squared distance, inertia and percentage of inertia
    ##########################################################################################################
    # Rows square distance to origin
    row_sqdisto = sum((Z .^2) .* transpose(col_weights),dims=2)

    # Rows inertia
    row_inertia = row_weights .* row_sqdisto

    # Rows percentage of inertia
    row_inertia_pct = 100*row_inertia /sum(row_inertia)

    # convert to dataframe
    row_infos = hcat(row_names,DataFrame([row_weights row_sqdisto row_inertia row_inertia_pct],["Weight","Sq. Dist","Inertia","% Inertia"]))
    
    ##########################################################################################
    ## Columns informations : weights, squared distance, inertia and percentage of inertia
    ##########################################################################################
    # Columns square distance to origin
    col_sqdisto = transpose(sum((Z .^2) .* row_weights,dims=1))

    # Columns inertia
    col_inertia = col_weights .* col_sqdisto

    # Columns percentage of inertia
    col_inertia_pct = 100*col_inertia/sum(col_inertia)

    #convert to dataframe
    col_infos = hcat(col_names,DataFrame([col_weights col_sqdisto col_inertia col_inertia_pct],["Weight","Sq. Dist","Inertia","% Inertia"]))

    # Generalized singular value decomposition (GSVD)
    svd_ = svd_triplet(Z;row_weights=row_weights,col_weights=col_weights,n_components=n_components)

    # Initialize NamedTuple
    res = (; :svd => svd_)

    # eigen values
    eigen_values = svd_.vs[1:max_components] .^2
    difference = vcat(reverse(diff(reverse(eigen_values))),missing)
    proportion = 100*eigen_values/sum(eigen_values)
    cumulative = cumsum(proportion)

    # Convert to dataframe
    eig_ = DataFrame((; :Dimensions => ["Dim."*string(i) for i in 1:max_components],
                        :Eigenvalue => eigen_values,
                        :Difference => difference, 
                        :Proportion => proportion, 
                        Symbol("Cum. Proportion") => cumulative))
    
    # Update namedtuple
    res = @insert res.eig = eig_;
    
    #######################################################################################################
    ## Rows informations : factor coordinates, contributions and square cosinus
    #######################################################################################################
    # Rows factor coordinates
    row_coord = svd_.U * LinearAlgebra.diagm(sqrt.(eigen_values)[1:n_components])
    
    # Rows contributions
    row_contrib = 100*((row_coord .^2) .* row_weights)./transpose(eigen_values[1:n_components])

    # Rows square cosinus
    row_cos2 = (row_coord .^2) ./ row_sqdisto

    # Store all informations
    row_ = (; :coord => hcat(row_names,DataFrame(row_coord,dim_index)),
              :cos2 => hcat(row_names,DataFrame(row_cos2,dim_index)),
              :contrib => hcat(row_names,DataFrame(row_contrib,dim_index)),
              :infos => row_infos)

    # Update namedtuple
    res = @insert res.row = row_;

    #####################################################################################################
    ## Columns informations : factor coordinates, contributions and square cosinus
    #####################################################################################################
    # Columns factor coordinates
    col_coord = svd_.V * LinearAlgebra.diagm(sqrt.(eigen_values)[1:n_components])

    # Columns contributions
    col_contrib = 100*((col_coord .^2) .* col_weights)./transpose(eigen_values[1:n_components])

    # Columns square cosinus
    col_cos2 = (col_coord .^2) ./ col_sqdisto

    # store all informations
    col_ = (; :coord => hcat(col_names,DataFrame(col_coord,dim_index)),
              :cos2 => hcat(col_names,DataFrame(col_cos2,dim_index)),
              :contrib => hcat(col_names,DataFrame(col_contrib,dim_index)),
              :infos => col_infos)

    # Update namedtuple
    res = @insert res.col = col_;
    return res
end;
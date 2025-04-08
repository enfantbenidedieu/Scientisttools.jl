"""
    CA(X;...)

Correspondence Analysis (CA)

...
# Description
Performs Correspondence Analysis (CA) over the data given in DataFrame X including supplementary row and/or column points, supplementary quantitative variables and supplementary categorical variables.

# Required arguments
- `X`: DataFrame of shape (n_rows, n_columns)
    Training data, where `n_rows` is the number of rows and `n_columns` is the number of columns.
    `X` is a contingency table containing absolute frequencies.

# Optional (default) arguments
- `n_components`: an integer specified the number of dimensions kept in the results (by default 5).
- `row_weights`: an optional rows weights (by default, a vector of 1 for uniform rows weights), the weights are given only for active rows.
- `row_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary rows.
- `col_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary columns.
- `quanti_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary quantitative variables.
- `quali_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary qualitative variables.
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is used as index (rows) names
    - If `false`: index (rows) names are created.

# Returns
A NamedTuple with following keys:
- `model`: string specifying the model fitted = 'ca'.
- `call`: namedtuple with some informations.
- `eig`: dataFrame containing all the eigenvalues, the difference between each eigenvalue, the percentage of variance and the cumulative percentage of variance.
- `svd`: namedtuple with results for generalied singular value decomposition (GSVD)
- `row`: namedtuple of DataFrame containing all the results for the active rows (factor coordinates, square cosinus, relative contributions, additionals informations).
- `col`: namedtuple of DataFrame containing all the results for the active columns (factor coordinates, square cosinus, relative contributions, additionals informations).
- `others`: others informations for Correspondence Analysis (kaiser threshold, etc.).
- `row_sup`: namedtuple of DataFrame containing all the results for the supplementary rows (factor coordinates, square cosinus, square distance to origin).
- `col_sup`: namedtuple of DataFrame containing all the results for the supplementary columns (factor coordinates, square cosinus, square distance to origin).
- `quanti_sup`: namedtuple of DataFrame containing all the results for the supplementative quantitative variables (factor coordinates, correlation between variables and axes, square cosinus).
- `quali_sup`: namedtuple of DataFrame containing all the results for supplementary qualitative variables (factor coordinates of each categories of each variables, vtest which is a criterion with a Normal distribution, and eta2 which is the square correlation coefficient between qualitative variable and an axe).
- `summary_quanti_sup`: dataframe specifying summary statistics for supplementary quantitative variables.
- `summary_quali_sup`: dataframe specifying summary statistics for supplementary qualitative variables.

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

# References
- Escofier B, Pagès J (2023), Analyses Factorielles Simples et Multiples. 5ed, Dunod
- Husson, F., Le, S. and Pages, J. (2009). Analyse de donnees avec R, Presses Universitaires de Rennes.
- Husson, F., Le, S. and Pages, J. (2010). Exploratory Multivariate Analysis by Example Using R, Chapman and Hall.
- Lebart L., Piron M., & Morineau A. (2006). Statistique exploratoire multidimensionnelle. Dunod, Paris 4ed.
- Pagès J. (2013). Analyse factorielle multiple avec R : Pratique R. EDP sciences
- Rakotomalala R. (2020), Pratique des méthodes factorielles avec Python, Université Lumière Lyon 2, Version 1.0

# See Also
`get_ca_row`, `get_ca_col`, `get_ca`, `dimdesc`, `predictCA`, `supvarCA`, `fviz_ca_row`, `fviz_ca_col`, `fviz_ca_biplot`.

...
"""
function CA(X::AbstractDataFrame;
            n_components::Int=5,
            row_weights=nothing,
            row_sup=nothing,
            col_sup=nothing,
            quanti_sup=nothing,
            quali_sup=nothing,
            first_col_as_index::Bool=true)::NamedTuple

    ##########################################################################################
    ## Check if X is a DataFrame
    ##########################################################################################
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a dataframe.")) end;
    
    ##########################################################################################
    ## Check if supplementary rows
    ##########################################################################################
    if row_sup !== nothing
        if isa(row_sup,Int)
            row_sup_idx = [Int(row_sup)]
        elseif isa(row_sup,AbstractVector)
            row_sup_idx = reduce(vcat, row_sup)
        elseif isa(row_sup,AbstractUnitRange)
            row_sup_idx = reduce(vcat,collect.([row_sup]))
        end;
    else
        row_sup_idx = nothing
    end;

    ##########################################################################################
    ## Check if supplementary columns
    ##########################################################################################
    if col_sup !== nothing
        if isa(col_sup,Int)
            col_sup_idx = [Int(col_sup)]
        elseif isa(col_sup,AbstractVector)
            col_sup_idx = reduce(vcat, col_sup)
        elseif isa(col_sup,AbstractUnitRange)
            col_sup_idx = reduce(vcat,collect.([col_sup]))
        end;
        col_sup_label = names(X)[col_sup_idx]
    else
        col_sup_label = nothing
    end;
    
    ##########################################################################################
    ## Check if supplementary quantitatives variables
    ##########################################################################################
    if quanti_sup !== nothing
        if isa(quanti_sup,Int)
            quanti_sup_idx = [Int(quanti_sup)]
        elseif isa(quanti_sup,AbstractVector)
            quanti_sup_idx = reduce(vcat, quanti_sup)
        elseif isa(quanti_sup,AbstractUnitRange)
            quanti_sup_idx = reduce(vcat,collect.([quanti_sup]))
        end;
        quanti_sup_label = names(X)[quanti_sup_idx]
    else
        quanti_sup_label = nothing
    end;
    
    ##########################################################################################
    ## Check if supplementary qualitatives variables
    ##########################################################################################
    if quali_sup !== nothing
        if isa(quali_sup,Int)
            quali_sup_idx = [Int(quali_sup)]
        elseif isa(quali_sup,AbstractVector)
            quali_sup_idx = reduce(vcat, quali_sup)
        elseif isa(quali_sup,AbstractUnitRange)
            quali_sup_idx = reduce(vcat,collect.([quali_sup]))
        end;
        quali_sup_label = names(X)[quali_sup_idx]
    else
        quali_sup_label = nothing
    end;

    ##########################################################################################
    ## Store a copy
    ##########################################################################################
    Xtot = copy(X)

    ##########################################################################################
    ## Drop supplementary qualitatives variables
    ##########################################################################################
    if quali_sup !== nothing
        X = X[!,Not(quali_sup_label)]
    end;

    ##########################################################################################
    ## Drop supplementary quantitatives variables
    ##########################################################################################
    if quanti_sup !== nothing
        X = X[!,Not(quanti_sup_label)]
    end;

    ##########################################################################################
    ## Drop supplementary columns
    ##########################################################################################
    if col_sup !== nothing
        X = X[!,Not(col_sup_label)]
    end;

    ##########################################################################################
    ## Drop supplementary rows
    ##########################################################################################
    if row_sup !== nothing
        # Extract supplementary rows
        X_row_sup = X[row_sup_idx,:]
        X = X[Not(row_sup_idx), :]
    end;

    ##########################################################################################
    ## Set active row labels
    ##########################################################################################
    if first_col_as_index
        row_names = DataFrame(Modalities=string.(X[!,1]))
        # Drop first columns
        X = X[!,2:end]
    else
        row_names = DataFrame(Modalities=["row" * string(i) for i in 1:nrow(X)])
    end;

    ##########################################################################################
    ## Create label for active columns
    ##########################################################################################
    col_names = DataFrame(Modalities=string.(names(X))) 

    # size of active elements
    n_rows, n_cols = nrow(row_names), nrow(col_names)

    # Convert data to matrix
    vquant = Matrix{Int}(X);

    ##########################################################################################
    ## Set row weights
    ##########################################################################################
    if row_weights === nothing
        row_weights = ones(n_rows)
    elseif !isa(row_weights,Array{<:Number,1})
        throw(ArgumentError("row_weights must be a vector of number."))
    elseif length(row_weights) != n_rows
        throw(DimensionMismatch("row_weights must be a vector with length $n_rows."))
    end;

    # Weighted X with the row weight
    wX = vquant .* row_weights

    # Apply row weights and sum
    total = sum(wX)

    # Table of frequencies
    freq =  wX /total

    # Compute row marge and columns marges
    col_marge, row_marge = vec(sum(freq,dims=1)), vec(sum(freq,dims=2))

    # Creation of matrix Z (useful to SVD)
    Z = ((freq ./ row_marge) ./ transpose(col_marge)) .- 1
    
    # QR decomposition (to set maximum number of components)
    qr_dec = qr(Matrix(Z))
    max_components = min(min(rank(Matrix(qr_dec.Q)), rank(Matrix(qr_dec.R))),n_rows - 1, n_cols - 1)
    
    # set number of components
    if n_components === nothing
        n_components = max_components
    elseif !isa(n_components,Int)
        throw(ArgumentError("n_components must be an integer."))
    elseif n_components < 1
        throw(ArgumentError("n_components must be equal or greater than 1."))
    else
        n_components = min(n_components,max_components)
    end;

    # set factor coordinates columns
    dim_index = ["Dim."*string(i) for i in 1:n_components]

    # store call informations
    call_ = (; :Xtot => Xtot, 
               :X => hcat(row_names,X), 
               :Z => hcat(row_names,DataFrame(Z,col_names[!,1])),
               :row_weights => row_weights, 
               :row_marge => row_marge,
               :col_marge => col_marge, 
               :n_components => n_components,
               :dim_index => dim_index,
               :row_sup => row_sup_idx,
               :col_sup => col_sup_label,
               :quanti_sup => quanti_sup_label,
               :quali_sup => quali_sup_label,
               :first_col_as_index => first_col_as_index)

    # Initialize named tuple
    res = (; :model => "ca", :call => call_);

    ##########################################################################################
    ## Fit Factor Analysis model
    ##########################################################################################
    fit_ = fitFA(Z,max_components,n_components=n_components,row_weights=row_marge,row_names=row_names,col_weights=col_marge,col_names=col_names);

    # Extract elements and update namedtuple
    res = @insert res.svd = fit_.svd;
    res = @insert res.eig = fit_.eig;

    # Update row and columns informations
    row_infos, col_infos = fit_.row.infos, fit_.col.infos
    rename!(row_infos,:Weight => :Margin)
    rename!(col_infos,:Weight => :Margin)
    insertcols!(row_infos, 2, :Weight => row_weights)

    # Set 
    row_ = (; :coord => fit_.row.coord, :contrib => fit_.row.contrib, :cos2 => fit_.row.cos2, :infos => row_infos)
    col_ = (; :coord => fit_.col.coord, :contrib => fit_.col.contrib, :cos2 => fit_.col.cos2, :infos => col_infos)
    
    # Update NamedTuple
    res = @insert res.row = row_;
    res = @insert res.col = col_;

    ###########################################################################################
    # Compute others indicators
    ###########################################################################################
    # Compute chi - squared test
    statistic,p_value,dof,expected_freq = chi2_contingency(wX)

    # Absolute residuals
    resid = wX - expected_freq

    # Standardized resid
    standardized_resid = resid ./ sqrt.(expected_freq)

    # Adjusted residuals
    adjusted_resid = (standardized_resid ./ sqrt.(1 .- row_marge)) ./transpose(sqrt.(1 .- col_marge))

    # Chi2 contribution
    chi2_contrib = 100*(standardized_resid .^2) ./ statistic

    # Attraction - Repulsion index
    iar = wX ./ expected_freq

    # Chi2 test
    chi2_test = DataFrame(Dict(:Names => ["statistic","dof","pvalue"], :value => [statistic,dof,p_value]))

    # Log-likelihood test
    gtest = g_test(wX)
    gtest = DataFrame(Dict(:Names => ["statistic","dof","pvalue"], :value => [gtest.statistic,gtest.dof,gtest.pvalue]))

    # Association test
    association = association_measure(vquant)

    # Kaiser threshold
    kaiser_threshold = mean(fit_.eig[!,2])
    kaiser_proportion_threshold = 100/max_components

    # Store others informations
    others_ = (; :resid => hcat(row_names,DataFrame(resid,col_names[!,1])),
                 :chi2 => chi2_test,
                 :g_test => gtest,
                 :resid => hcat(row_names,DataFrame(resid,col_names[!,1])),
                 :resid_adj => hcat(row_names,DataFrame(adjusted_resid,col_names[!,1])),
                 :rstandard => hcat(row_names,DataFrame(standardized_resid,col_names[!,1])),
                 :chi2_contrib => hcat(row_names,DataFrame(chi2_contrib,col_names[!,1])),
                 :attraction => hcat(row_names,DataFrame(iar,col_names[!,1])),
                 :association => association,
                 :kaiser => (; :threshold => kaiser_threshold, :proportion_threshold => kaiser_proportion_threshold))

    # Update namedtuple
    res = @insert res.others = others_;

    ##############################################################################
    ## Compute statistics for supplementary rows
    ##############################################################################
    if row_sup !== nothing
        # Set supplementary rows names
        if first_col_as_index
            row_sup_names = DataFrame(Modalities=string.(X_row_sup[!,1]))
            # Drop first columns
            X_row_sup = X_row_sup[!,2:end] 
        else
            row_sup_names = DataFrame(Modalities=["row" * string(i) for i in (n_rows+1):(n_rows+nrow(X_row_sup))])
        end;

        # Convert to matrice
        X_row_sup = Matrix{Int}(X_row_sup)

        # Rows sum
        row_sum = vec(sum(X_row_sup,dims=2))

        # Standardize with the row sum
        Z_row_sup = X_row_sup ./ row_sum

        # Supplementary rows factor coordinates
        row_sup_coord = Z_row_sup * fit_.svd.V[:,1:n_components]

        # Supplementary rows square distance to origin
        row_sup_sqdisto = sum(((Z_row_sup .- transpose(col_marge)) .^2) ./ transpose(col_marge),dims=2)

        # Supplementary rows cos2
        row_sup_cos2 = (row_sup_coord .^2) ./ row_sup_sqdisto

        # Store all informations
        row_sup_ = (; :coord => hcat(row_sup_names,DataFrame(row_sup_coord,dim_index)),
                      :cos2 => hcat(row_sup_names,DataFrame(row_sup_cos2,dim_index)),
                      :dist2 => hcat(row_sup_names,DataFrame("Sq. Dist" => vec(row_sup_sqdisto))))

        # Update namedtuple
        res = @insert res.row_sup = row_sup_;
    end;

    ################################################################################################
    ## Compute statistics for supplementary columns
    ################################################################################################
    if col_sup !== nothing
        # Select supplementary columns
        X_col_sup = Xtot[!,col_sup_label]

        # Transform to DataFrame if Vector
        if isa(X_col_sup,AbstractVector) X_col_sup = DataFrame(col_sup_label => X_col_sup) end;

        # Remove supplementary rows
        if row_sup !== nothing X_col_sup = X_col_sup[Not(row_sup_idx),:] end;

        # Set supplementary columns names
        col_sup_names = DataFrame(Modalities=string.(col_sup_label))

        # Weighted with rows weights
        X_col_sup = Matrix{Int}(X_col_sup) .* row_weights

        # Columns sum
        col_sum = vec(sum(X_col_sup,dims=1))

        # Standardize with columns sum
        Z_col_sup = X_col_sup ./ transpose(col_sum)

        # Supplementary columns factor coordinates
        col_sup_coord = transpose(Z_col_sup) * fit_.svd.U[:,1:n_components]

        # Supplementary columnns square distance to origin
        col_sup_sqdisto = vec(sum(((Z_col_sup .- row_marge) .^2)./row_marge,dims=1))

        # Supplementary columns cos2
        col_sup_cos2 = (col_sup_coord .^2) ./ col_sup_sqdisto

        # Store all informations
        col_sup_ = (; :coord => hcat(col_sup_names,DataFrame(col_sup_coord,dim_index)),
                      :cos2 => hcat(col_sup_names,DataFrame(col_sup_cos2,dim_index)),
                      :dist2 => hcat(col_sup_names,DataFrame("Sq. Dist" => vec(col_sup_sqdisto))))

        # Update namedtuple
        res = @insert res.col_sup = col_sup_;
    end;

    #######################################################################################
    ## Compute statistics for supplementary quantitative variables
    #######################################################################################
    if quanti_sup !== nothing
        # Select supplementary quantitatives variables
        X_quanti_sup = Xtot[!,quanti_sup_label]

        # Transform to DataFrame if Vector
        if isa(X_quanti_sup,AbstractVector) X_quanti_sup = DataFrame(quanti_sup_label => X_quanti_sup) end;
        
        # Remove supplementary rows
        if row_sup !== nothing X_quanti_sup = X_quanti_sup[Not(row_sup_idx),:] end;

        # Store all informations
        quanti_sup_ = predict_quanti_sup(X_quanti_sup,row_marge,fit_.svd.U)

        # Update namedtuple
        res = @insert res.quanti_sup = quanti_sup_.quanti_sup;
        res = @insert res.summary_quanti_sup = quanti_sup_.summary;
    end;

    ##########################################################################################
    ## Compute statistics for supplementary qualitative variables
    ##########################################################################################
    if quali_sup !== nothing
        # Select supplementary qualitatives variables
        X_quali_sup = Xtot[!,quali_sup_label]
        
        # Transform to DataFrame if Vector
        if isa(X_quali_sup,AbstractVector) X_quali_sup = DataFrame(quali_sup_label => X_quali_sup) end;
        
        # Remove supplementary rows
        if row_sup !== nothing X_quali_sup = X_quali_sup[Not(row_sup_idx),:] end;
        
        # Statistics for qualitatives variables
        summary_quali_sup = freq_table(X_quali_sup)

        # Set supplementary categories names
        mod_sup_names = select(summary_quali_sup,:Modalities)

        # Compute data
        quali_sup_data = reduce(vcat,[sum_table(X,X_quali_sup,col) for col in names(X_quali_sup)])
        
        # Convert to matrice
        quali_sup_data = Matrix{Int}(quali_sup_data[!,2:end])
        
        # Standardize with the row sum
        Z_quali_sup = quali_sup_data ./ vec(sum(quali_sup_data,dims=2))

        # Supplementary categories factor coordinates
        quali_sup_coord = Z_quali_sup * fit_.svd.V[:,1:n_components]
        
        # Supplementary categories square distance to origin
        quali_sup_sqdisto = sum(((Z_quali_sup .- transpose(col_marge)) .^2) ./ transpose(col_marge),dims=2)

        # Supplementary categories cos2
        quali_sup_cos2 = (quali_sup_coord .^2) ./ quali_sup_sqdisto

        # Create disjonctif table
        dummies = Matrix{Int}(get_dummies(X_quali_sup))
        # Apply correction : weighted count by categories
        n_k = vec(sum(dummies .* row_marge,dims=1) .* total)

        # Supplementary categories value-test
        if total > 1
            quali_sup_vtest = quali_sup_coord .* sqrt.((n_k .* (total - 1)) ./ (total .- n_k))
        else
            quali_sup_vtest = quali_sup_coord .* sqrt.(n_k)
        end;

        # Supplementary qualitative variables square correlation ratio
        quali_sup_eta2 = reduce(vcat,[function_eta2(X_quali_sup,lab,fit_.row.coord[!,2:end];w=row_marge) for lab in names(X_quali_sup)])
        
        # Store all informations
        quali_sup_ = (; :coord => hcat(mod_sup_names,DataFrame(quali_sup_coord,dim_index)),
                        :cos2 => hcat(mod_sup_names,DataFrame(quali_sup_cos2,dim_index)),
                        :vtest => hcat(mod_sup_names,DataFrame(quali_sup_vtest,dim_index)),
                        :eta2 => quali_sup_eta2,
                        :dist2 => hcat(mod_sup_names,DataFrame("Sq. Dist" => vec(quali_sup_sqdisto))))

        # Update namedtuple
        res = @insert res.quali_sup = quali_sup_;
        res = @insert res.summary_quali_sup = summary_quali_sup;
    end;

    return res
end;

"""
    predictCA(self, X;...)

Predict projection for new rows with Correspondence Analysis (CA)

....
# Description
Performs the coordinates, square cosinus and square distance to origin of new rows with Correspondence Analysis (CA).

# Required arguments
- `self`: a CA NamedTuple
- `X`: a DataFrame in which to look for columns with which to predict. X must contain columns with the same names as the original data.

# Optionals (default) arguments
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is used as index (rows) names;
    - If `false`: index (rows) names are created.

# Returns
A NamedTuple of DataFrame containing all the results for the new rows including:
- `coord`: factor coordinates (scores) of the new rows;
- `cos2`: square cosinus of the new rows;
- `dist2`: square distance to origin of the new rows.

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Predict supplementary rows
julia> donnee = get_dataset("children",choice="row_sup");
julia> row_sup = predictCA(res_ca,donnee);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function predictCA(self::NamedTuple,
                   X::AbstractDataFrame;
                   first_col_as_index::Bool=true)::NamedTuple

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;

    # Check if X is a DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;
    
    # Check if CA model
    if self.model != "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;

    # Set row names
    if first_col_as_index
        row_names = DataFrame(Modalities=string.(X[!,1]))
        # Drop first columns
        X = X[!,2:end]
    else
        row_names = DataFrame(Modalities=["row" * string(i) for i in (nrow(self.call.X)+1):(nrow(self.call.X)+nrow(X))])
    end;

    # Check if columns are aligned
    if ncol(X) != ncol(self.call.X[!,2:end]) throw(ArgumentError("Columns aren't aligned.")) end;

    # Convert to matrice
    X = Matrix{Int}(X)

    # Rows sum
    row_sum = vec(sum(X,dims=2))

    # Standardize with the row sum
    Z = X ./ row_sum

    # Factor coordinates
    coord = Z * self.svd.V[:,1:self.call.n_components]

    # Square distance to origin
    sqdisto = sum(((Z .- transpose(self.call.col_marge)) .^2) ./ transpose(self.call.col_marge),dims=2)

    # Square cosinus
    cos2 = (coord .^2) ./ sqdisto

    # Store all informations
    res = (; :coord => hcat(row_names,DataFrame(coord,self.call.dim_index)),
             :cos2 => hcat(row_names,DataFrame(cos2,self.call.dim_index)),
             :dist2 => hcat(row_names,DataFrame("Sq. Dist" => vec(sqdisto))));
    return res
end;

"""
    supvarCA(self;...)

Supplementary columns/variables with Correspondence Analysis (CA)

...
# Description
Performns the coordinates, square cosinus and square distance to origin of supplementary columns/variables with Correspondence Analysis (CA).

# Required arguments
- `self`: a CA NamedTuple

# Optional (default) arguments
- `X_col_sup`: a DataFrame of supplementary columns (default = nothing)
- `X_quanti_sup`: a DataFrame of supplementary quantitative variables (default = nothing)
- `X_quali_sup`: a DataFrame of supplementary qualitative variables (default = nothing)
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is removed from DataFrame;
    - If `false`: the first columns is kepted in the DataFrame.

# Returns
A NamedTuple of NamedTuple containing all the results for the new variables including:
- `col_sup`: NamedTuple of DataFrame containing all the results for the new columns including:
    - `coord`: factor coordinates (scores) of the new columns;
    - `cos2`: square cosinus of the new columns;
    - `dist2`: square distance to origin of the new columns.

- `quanti_sup`: NamedTuple of DataFrame containing the results for the supplementary quantitative variables including:
    - `coord`: factor coordinates (scores) for the supplementary quantitative variables;
    - `cos2`: square cosinus for the supplementary quantitative variables;

- `quali_sup`: NamedTuple of DataFrame containing the results for the supplementary qualitative variables including:
    - `coord`: factor coordinates (scores) for the supplementary categories
    - `cos2`: square cosinus for the supplementary categories
    - `vtest`: value-test for the supplementary categories
    - `dist2`: square distance to origin for the supplementary categories
    - `eta2`: square correlation ratio for the supplementary qualitative variables.

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Predict supplementary variables
julia> col_sup = get_dataset("children",choice="col_sup");
julia> quali_sup = get_dataset("children",choice="quali_sup");
julia> sup_var = supvarCA(res_pca,X_col_sup=col_sup,X_quanti_sup=col_sup,X_quali_sup=quali_sup);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function supvarCA(self::NamedTuple;
                  X_col_sup=nothing,
                  X_quanti_sup=nothing,
                  X_quali_sup=nothing,
                  first_col_as_index::Bool=true)::NamedTuple

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    
    # Check if CA model
    if self.model != "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;
        
    # Check if all are empty
    if X_col_sup === nothing && X_quanti_sup === nothing && X_quali_sup === nothing  throw(ArgumentError("At least one shouldn't be empty."))  end;

    ##############################################################################################
    ## Statistics for supplementary columns
    ##############################################################################################
    if X_col_sup !== nothing
        # Check if X_col_sup is a DataFrame
        if !isa(X_col_sup,AbstractDataFrame) throw(ArgumentError("X_col_sup must be a DataFrame.")) end;
        
        # Check if first columns is a vector of string
        if first_col_as_index X_col_sup = X_col_sup[!,2:end] end;

        # Set new columns names
        col_sup_names = DataFrame(Modalities=string.(names(X_col_sup)))

        # Weighted with rows weights
        X_col_sup = Matrix{Int}(X_col_sup) .* self.call.row_weights

        # Columns sum
        col_sum = vec(sum(X_col_sup,dims=1))

        # Standardize with columns sum
        Z_col_sup = X_col_sup ./ transpose(col_sum)

        # Supplementary columns factor coordinates
        col_sup_coord = transpose(Z_col_sup) * self.svd.U[:,1:self.call.n_components]

        # Supplementary columnns square distance to origin
        col_sup_sqdisto = vec(sum(((Z_col_sup .- self.call.row_marge) .^2)./self.call.row_marge,dims=1))

        # Square cosinus
        col_sup_cos2 = (col_sup_coord .^2) ./ col_sup_sqdisto

        # Store all informations
        col_sup_ = (; :coord => hcat(col_sup_names,DataFrame(col_sup_coord,self.call.dim_index)),
                      :cos2 => hcat(col_sup_names,DataFrame(col_sup_cos2,self.call.dim_index)),
                      :dist2 => hcat(col_sup_names,DataFrame("Sq. Dist" => vec(col_sup_sqdisto))));
    else
        col_sup_ = nothing;
    end;

    #######################################################################################################
    ## Statistics for supplementary quantitative variables
    #######################################################################################################
    if X_quanti_sup !== nothing
        # Check if X_quanti_sup is a DataFrame
        if !isa(X_quanti_sup,AbstractDataFrame) throw(ArgumentError("X_quanti_sup must be a DataFrame.")) end;
        
        # Check if first columns is a vector of string
        if first_col_as_index X_quanti_sup = X_quanti_sup[!,2:end] end;

        # Store all informations
        quanti_sup_ = predict_quanti_sup(X_quanti_sup,self.call.row_marge,self.svd.U);
    else
        quanti_sup_ = nothing;
    end;

    ################################################################################################################
    ## Statistics for supplementary qualitative variables
    ################################################################################################################
    if X_quali_sup !== nothing
        # Check if X_quali_sup is a DataFrame
        if !isa(X_quali_sup,AbstractDataFrame) throw(ArgumentError("X_quali_sup must be a DataFrame.")) end;
        
        # Check if first columns is a vector of string
        if first_col_as_index X_quali_sup = X_quali_sup[!,2:end] end;

        # Compute data
        quali_sup_data = reduce(vcat,[sum_table(self.call.X[!,2:end],X_quali_sup,col) for col in names(X_quali_sup)])

        # Set supplementary categories names
        mod_sup_names = select(quali_sup_data,:Modalities)
        
        # Convert to matrice
        quali_sup_data = Matrix{Int}(quali_sup_data[!,2:end])
        
        # Standardize with the row sum
        Z_quali_sup = quali_sup_data ./ vec(sum(quali_sup_data,dims=2))

        # Fcator coordinates
        quali_sup_coord = Z_quali_sup * self.svd.V[:,1:self.call.n_components]
        
        # Square distance to origin
        quali_sup_sqdisto = sum(((Z_quali_sup .- transpose(self.call.col_marge)) .^2) ./ transpose(self.call.col_marge),dims=2)

        # Square cosinus
        quali_sup_cos2 = (quali_sup_coord .^2) ./ quali_sup_sqdisto

        # sum of elements
        total = sum(quali_sup_data)

        # Create disjonctif table
        dummies = Matrix{Int}(get_dummies(X_quali_sup))
        # Apply correction : weighted count by categories
        n_k = vec(sum(dummies .* self.call.row_marge,dims=1) .* total)

        # Supplementary categories value-test
        if total > 1
            quali_sup_vtest = quali_sup_coord .* sqrt.((n_k .* (total - 1)) ./ (total .- n_k))
        else
            quali_sup_vtest = quali_sup_coord .* sqrt.(n_k)
        end;

        # Supplementary qualitative variables square correlation ratio
        quali_sup_eta2 = reduce(vcat,[function_eta2(X_quali_sup,lab,self.row.coord[!,2:end],w=self.call.row_marge) for lab in names(X_quali_sup)])
        
        # Store all informations
        quali_sup_ = (; :coord => hcat(mod_sup_names,DataFrame(quali_sup_coord,self.call.dim_index)),
                        :cos2 => hcat(mod_sup_names,DataFrame(quali_sup_cos2,self.call.dim_index)),
                        :vtest => hcat(mod_sup_names,DataFrame(quali_sup_vtest,self.call.dim_index)),
                        :eta2 => quali_sup_eta2,
                        :dist2 => hcat(mod_sup_names,DataFrame("Sq. Dist" => vec(quali_sup_sqdisto))));
    else
        quali_sup_ = nothing;
    end;

    return (; :col_sup => col_sup_, :quanti_sup => quanti_sup_, :quali_sup => quali_sup_)
end;
# Principal Component Analysis

#=
"""
    PCA(X;...)

Principal Components Analysis (PCA)

...
# Description

Performs Principal Component Analysis (PCA) over the data given in DataFrame X with supplementary individuals, supplementary quantitative variables and supplementary categorical variables. 

Missing values are replaced by the column mean(default), by median or dopped.

# Required arguments
- `X`: DataFrame of shape (`n_samples`, `n_columns`)
    Training data, where `n_samples` is the number of samples and `n_columns` is the number of columns.

# Optional (default) arguments
- `standardize`: a boolean, default = true:
    - If `true`: the data are scaled to unit variance;
    - If `false`: the data are not scaled to unit variance.
- `n_components`: an integer specified the number of dimensions kept in the results (by default 5).
- `ind_weights`: an optional individuals weights (by default, a vector of 1/(number of active individuals) for uniform individuals weights), the weights are given only for active individuals.
- `var_weights`: an optional variables weights (by default, a vector of 1 for uniform variables weights), the weights are given only for active variables.
- `ind_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary individuals.
- `quanti_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary quantitative variables.
- `quali_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary qualitative variables.
- `fill_missing`: a symbol specified method for missing treatment. Allowed values are:
    - `:mean` for fill missing with mean
    - `:median` for fill missing with median
    - `drop` for dropping missing values.
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is used as index (individuals) names
    - If `false`: index (individuals) names are created.

# Returns
A NamedTuple with following keys:
- `model`: string specifying the model fitted = 'pca'.
- `call`: namedtuple with some informations.
- `eig`: dataFrame containing all the eigenvalues, the difference between each eigenvalue, the percentage of variance and the cumulative percentage of variance.
- `svd`: namedtuple with results for generalied singular value decomposition (GSVD)
- `ind`: namedtuple of DataFrame containing all the results for the active individuals (factor coordinates, square cosinus, relative contributions, additionals informations).
- `var`: namedtuple of DataFrame containing all the results for the active variables (factor coordinates, square cosinus, relative contributions, additionals informations).
- `others`: others informations for Principal Components Analysis (kaiser threshold, etc.).
- `ind_sup`: namedtuple of DataFrame containing all the results for the supplementary individuals (factor coordinates, square cosinus,square distance to origin).
- `quanti_sup`: namedtuple of DataFrame containing all the results for the supplementative quantitative variables (factor coordinates, correlation between variables and axes, square cosinus).
- `quali_sup`: namedtuple of DataFrame containing all the results for supplementary qualitative variables (factor coordinates of each categories of each variables, vtest which is a criterion with a Normal distribution, and eta2 which is the square correlation coefficient between qualitative variable and an axe).
- `summary_quanti`: dataframe specifying summary statistics for quantitative variables (actives and/or supplementary).
- `summary_quali`: dataframe specifying summary statistics for supplementary qualitative variables.
- `association`: dataframe specifying association test. If supplementary qualitative are greater than 2.
- `chi2_test`: dataframe specifying chi2 contingency test. If supplementary qualitative are greater than 2.

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

# References
- Bry X. (1996), Analyses factorielles multiple, Economica
- Bry X. (1999), Analyses factorielles simples, Economica
- Escofier B., Pagès J. (2023), Analyses Factorielles Simples et Multiples. 5ed, Dunod
- Saporta G. (2006). Probabilites, Analyse des données et Statistiques. Technip
- Husson, F., Le, S. and Pages, J. (2010). Exploratory Multivariate Analysis by Example Using R, Chapman and Hall.
- Lebart L., Piron M., & Morineau A. (2006). Statistique exploratoire multidimensionnelle. Dunod, Paris 4ed.
- Pagès J. (2013). Analyse factorielle multiple avec R : Pratique R. EDP sciences
- Rakotomalala, R. (2020). Pratique des méthodes factorielles avec Python. Université Lumière Lyon 2. Version 1.0
- Tenenhaus, M. (2006). Statistique : Méthodes pour décrire, expliquer et prévoir. Dunod.

# See Also
`get_pca_ind`, `get_pca_var`, `get_pca`, `dimdesc`, `reconst`, `predictPCA`, `supvarPCA`, `fviz_pca_ind`, `fviz_pca_var`, `fviz_pca_biplot`.

...
"""
=#
function PCA(X::AbstractDataFrame;
             standardize::Bool=true,
             n_components::Int=5,
             ind_weights=nothing,
             var_weights=nothing,
             ind_sup=nothing,
             quanti_sup=nothing,
             quali_sup=nothing,
             fill_missing::Symbol=:mean,
             first_col_as_index::Bool=true)::NamedTuple

    ###########################################################################################
    ## Check if X is a DataFrame
    ###########################################################################################
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    ##########################################################################################
    ## Check if supplementary individuals
    ##########################################################################################
    if ind_sup !== nothing
        if isa(ind_sup,Int)
            ind_sup_idx = [Int(ind_sup)]
        elseif isa(ind_sup,AbstractVector)
            ind_sup_idx = reduce(vcat, ind_sup)
        elseif isa(ind_sup,AbstractUnitRange)
            ind_sup_idx = reduce(vcat,collect.([ind_sup]))
        end;
    else
        ind_sup_idx = nothing
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

    ###########################################################################################
    ## Check if supplementary qualitatives variables
    ###########################################################################################
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

    #############################################################################################
    ## Fill missing with mean or median
    #############################################################################################
    if fill_missing == :drop
        X = dropmissing(X)
    elseif fill_missing in [:mean,:median]
        X = missing_imputation(X,method=fill_missing)
    else
        throw(ArgumentError("'fill_missing' should be one of 'drop', 'mean', 'median'"))
    end;

    #############################################################################################
    ## Make a copy of the data
    #############################################################################################
    Xtot = copy(X)

    #############################################################################################
    ## Drop supplementary qualitatives variables
    #############################################################################################
    if quali_sup !== nothing X = X[!,Not(quali_sup_label)] end;

    #############################################################################################
    ## Drop supplementary quantitatives variables
    #############################################################################################
    if quanti_sup !== nothing X = X[!,Not(quanti_sup_label)] end;
    
    #############################################################################################
    # Drop supplementary individuals
    #############################################################################################
    if ind_sup !== nothing
        # Extract supplementary individuals
        X_ind_sup = X[ind_sup_idx,:]
        X = X[Not(ind_sup_idx), :]
    end;

    #############################################################################################
    # Set active individuals names
    #############################################################################################
    if first_col_as_index
        ind_names = DataFrame(Individuals=string.(X[!,1]))
        # Drop first columns
        X = X[!,2:end]
    else
        ind_names = DataFrame(Individuals=["ind" * string(i) for i in 1:nrow(X)])
    end;

    ###############################################################################################
    # Set active variables names
    ###############################################################################################
    var_names = DataFrame(Variables=string.(names(X))) 

    # size of active elements
    n_rows, n_vars = nrow(ind_names), nrow(var_names)
    
    # Statistics of quantitatives variables
    summary_quanti = summary_stats(X)

    # Convert data to matrix
    vquant = Matrix{Float64}(X);

    #############################################################################################
    ## Set individuals weights
    #############################################################################################
    if ind_weights === nothing
        ind_weights = ones(n_rows)/n_rows
    elseif !isa(ind_weights,Array{<:Number,1})
        throw(ArgumentError("ind_weights must be a vector of number."))
    elseif length(ind_weights) != n_rows
        throw(DimensionMismatch("ind_weights must be a vector with length $n_rows."))
    else
        ind_weights = ind_weights/sum(ind_weights)
    end;

    #############################################################################################
    ## Set variables weights
    #############################################################################################
    if var_weights === nothing
        var_weights = ones(n_vars)
    elseif !isa(var_weights,Array{<:Number,1})
        throw(ArgumentError("var_weights must be a vector of number."))
    elseif length(var_weights) != n_vars
        throw(DimensionMismatch("var_weights must be a vector with length $n_vars."))
    else
        var_weights = vec(var_weights)
    end;

    # Average and standard deviation
    wmean, wstd = descrstatsw(vquant,weights(ind_weights));

    # set standardize
    if standardize
        stds = wstd
    else
        stds = reshape(ones(n_vars),(1,n_vars))
    end;

    # Standardisation : Z = (X - mu)/sigma
    Z = (vquant .- wmean) ./ stds

    # QR decomposition (to set maximum number of components)
    qr_dec = qr(Matrix(Z))
    max_components = min(min(rank(Matrix(qr_dec.Q)),rank(Matrix(qr_dec.R))), n_rows - 1, n_vars)
    
    ################################################################################################
    ## Set number of components
    ################################################################################################
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
    dim_index = ["Dim." * string(i) for i in 1:n_components]

    # store call informations
    call_ = (; :Xtot => Xtot, 
               :X => hcat(ind_names,X), 
               :Z => hcat(ind_names,DataFrame(Z,var_names[!,1])),
               :ind_weights => ind_weights, 
               :var_weights => var_weights,
               :means => wmean, 
               :std => stds, 
               :n_components => n_components,
               :dim_index => dim_index,
               :ind_sup => ind_sup_idx,
               :quanti_sup => quanti_sup_label,
               :quali_sup => quali_sup_label,
               :first_col_as_index => first_col_as_index)

    # Initialize named tuple
    res = (; :model => "pca", :call => call_);

    ##################################################################################################
    # Fit Factor Analysis model
    ##################################################################################################
    fit_ = fitFA(Z,max_components,n_components=n_components,row_weights=ind_weights,row_names=ind_names,col_weights=var_weights,col_names=var_names);

    # Extract elements and update namedtuple
    res = @insert res.svd = fit_.svd;
    res = @insert res.eig = fit_.eig;
    res = @insert res.ind = fit_.row;
    res = @insert res.var = fit_.col;

    ##################################################################################################
    ## Others metrics
    ##################################################################################################
    # Bartlett statistics
    bs_stats = - (n_rows - 1 - (2*n_vars + 5)/6)*sum(log.(fit_.eig[!,2]))
    bs_dof = n_vars*(n_vars - 1)/2 
    bs_pvalue = 1 - cdf(Chisq(bs_dof),bs_stats)
    bs_test = DataFrame(Dict(:Names => ["|CORR.MATRIX|","statistic","dof","p-value"], :value => [sum(log.(fit_.eig[!,2])),bs_stats,bs_dof,bs_pvalue]))
    
    # Kaiser threshold
    kaiser_threshold = mean(fit_.eig[!,2])
    kaiser_proportion_threshold = 100/max_components

    # KSS(Karlis - Saporta - Spinaki) threshold
    kss_threshold = 1 + 2*sqrt((n_vars - 1)/(n_rows - 1))

    # Broken stick threshold
    broken_stick_threshold = reverse(cumsum(reverse([1/i for i in 1:max_components])))
    broken_stick = hcat(select(fit_.eig,[:Dimensions,:Eigenvalue]),DataFrame("Threshold" => broken_stick_threshold))

    # Store others informations
    others_ = (; :bartlett => bs_test,
                 :kaiser => (; :threshold => kaiser_threshold, :proportion_threshold => kaiser_proportion_threshold),
                 :karlis_saporta_spinaki => kss_threshold,
                 :broken_stick_threshold => broken_stick)

    # Update namedtuple
    res = @insert res.others = others_;

    #######################################################################################################
    # Statistics for supplementary individuals
    #######################################################################################################
    if ind_sup !== nothing
        # Set supplementary individuals names
        if first_col_as_index
            ind_sup_names = DataFrame(Individuals=string.(X_ind_sup[!,1]))
            # Drop first columns
            X_ind_sup = X_ind_sup[!,2:end] 
        else
            ind_sup_names = DataFrame(Individuals=["ind" * string(i) for i in (n_rows+1):(n_rows+nrow(X_ind_sup))])
        end;

        # Standardisation
        Z_ind_sup = (Matrix{Float64}(X_ind_sup) .- wmean) ./ stds

        # Store all informations
        ind_sup_ = predict_ind_sup(Z_ind_sup,ind_sup_names,var_weights,fit_.svd.V)

        # Update namedtuple
        res = @insert res.ind_sup = ind_sup_;
    end;

    ######################################################################################################
    ## Statistics for supplementary quantitatives variables
    ######################################################################################################
    if quanti_sup !== nothing
        # Select supplementary quantitatives variables
        X_quanti_sup = Xtot[!,quanti_sup_idx]

        # Transform to DataFrame if Vector
        if isa(X_quanti_sup,AbstractVector) X_quanti_sup = DataFrame(quanti_sup_label => X_quanti_sup) end;

        # Remove supplementary individuals
        if ind_sup !== nothing X_quanti_sup = X_quanti_sup[Not(ind_sup_idx),:] end;

        # Store all informations
        quanti_sup_ = predict_quanti_sup(X_quanti_sup,ind_weights,fit_.svd.U)

        # Extract elements
        summary_quanti_sup = quanti_sup_.summary

        # Insert group columns
        insertcols!(summary_quanti,1, :group => "active")
        insertcols!(summary_quanti_sup,1, :group => "sup")

        # Concatenate
        summary_quanti = vcat(summary_quanti,summary_quanti_sup)
        
        # Update namedtuple
        res = @insert res.quanti_sup = quanti_sup_.quanti_sup;
    end;

    ########################################################################################
    ## Statistics for supplementary qualitatives variables
    ########################################################################################
    if quali_sup !== nothing
        # Select supplementary qualitatives variables
        X_quali_sup = Xtot[!,quali_sup_idx]

        # Transform to DataFrame if Vector
        if isa(X_quali_sup,AbstractVector) X_quali_sup = DataFrame(quali_sup_label => X_quali_sup) end;

        # Remove supplementary individuals
        if ind_sup !== nothing X_quali_sup = X_quali_sup[Not(ind_sup_idx),:] end;

        ## Compute association test
        if ncol(X_quali_sup) > 1
            # All combinations
            cols = collect(combinations(names(X_quali_sup),2))
            col_names = reduce(vcat,[DataFrame(Dict("variable1" => col[1],"variable2" => col[2])) for col in cols])
            cont_table = Dict(col => freqtable(X_quali_sup[!,col],Symbol(col[1]),Symbol(col[2])) for col in cols)

            # Chi2 statistic test
            chi2_test = hcat(col_names,reduce(vcat,[DataFrame([chi2_contingency(cont_table[col])[[:statistic,:dof,:pvalue]]]) for col in cols]))
            # Association test
            association = hcat(col_names,reduce(vcat,[association_measure(cont_table[col])[!,2:end] for col in cols]))

            # Update NamedTuple
            res = @insert res.chi2_test = chi2_test;
            res = @insert res.association = association;
        end;

        # Statistics for qualitatives variables
        summary_quali_sup = freq_table(X_quali_sup)

        # Set supplementary categories names
        mod_sup_names = select(summary_quali_sup,:Modalities)

        # Conditional weighted average with original dataset
        barycentre = conditional_weighted_average(X,X_quali_sup,w=ind_weights)

        # Standardization
        Z_quali_sup = (Matrix{Float64}(barycentre[!,2:end]) .- wmean) ./ stds

        # Supplementary categories factor coordinates
        quali_sup_coord = (Z_quali_sup .* transpose(var_weights)) * fit_.svd.V

        # Supplementary categories squared distance to origin
        quali_sup_sqdisto = sum((Z_quali_sup .^2) .* transpose(var_weights),dims=2)
        
        # Supplementary categories cos2
        quali_sup_cos2 = (quali_sup_coord .^2) ./quali_sup_sqdisto

        # Supplementary categories values-tests
        n_k = summary_quali_sup.Effectifs
        quali_sup_vtest = (quali_sup_coord ./transpose(sqrt.(fit_.eig[1:n_components,2]))) .* sqrt.(((n_rows - 1).*n_k)./(n_rows .- n_k));

        # Supplementary qualitative variables square correlation ratio
        quali_sup_eta2 = reduce(vcat,[function_eta2(X_quali_sup,lab,fit_.row.coord[!,2:end];w=ind_weights) for lab in names(X_quali_sup)])

        # Store all informations
        quali_sup_ = (; :barycentre => barycentre,
                        :coord => hcat(mod_sup_names,DataFrame(quali_sup_coord,dim_index)),
                        :cos2 => hcat(mod_sup_names,DataFrame(quali_sup_cos2,dim_index)),
                        :vtest => hcat(mod_sup_names,DataFrame(quali_sup_vtest,dim_index)),
                        :eta2 => quali_sup_eta2,
                        :dist2 => hcat(mod_sup_names,DataFrame("Sq. Dist" => vec(quali_sup_sqdisto))))

        # Update namedtuple
        res = @insert res.quali_sup = quali_sup_;
        res = @insert res.summary_quali_sup = summary_quali_sup;
    end;

    # Update with summary
    res = @insert res.summary_quanti = summary_quanti;

    return res
end;

"""
    predictPCA(self,X;...)

Predict projection for new individuals with Principal Components Analysis (PCA).

...
# Description
Performs the coordinates, square cosinus and square distance to origin of new individuals with Principal Component Analysis (PCA).

# Required arguments
- `self`: a PCA NamedTuple
- `X`: a DataFrame in which to look for variables with which to predict. X must contain columns with the same names as the original data.

# Optionals (default) arguments
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is used as index (individuals) names;
    - If `false`: index (individuals) names are created.

# Returns
A NamedTuple of DataFrame containing all the results for the new individuals including:
- `coord`: factor coordinates (scores) of the new individuals;
- `cos2`: square cosinus of the new individuals;
- `dist2`: square distance to origin of the new individuals.

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Predict supplementary individuals
julia> donnee = get_dataset("decathlon",choice="ind_sup");
julia> ind_sup = predictPCA(res_pca,donnee);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function predictPCA(self::NamedTuple,
                    X::AbstractDataFrame;
                    first_col_as_index::Bool=true)::NamedTuple

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;

    # Check if X is a DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Check if PCA model
    if self.model != "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;

    # Set individuals names
    if first_col_as_index 
        ind_names = DataFrame(Individuals=string.(X[!,1]))
        # Drop first columns
        X = X[!,2:end]
    else
        ind_names = DataFrame(Individuals=["ind" * string(i) for i in (nrow(self.call.X)+1):(nrow(self.call.X)+nrow(X))])
    end;

    # Check if columns are aligned
    if ncol(X) != ncol(self.call.X[!,2:end]) throw(ArgumentError("Columns aren't aligned.")) end;

    # Convert to float
    X = convert_to_float(X)

    # Creation of matrix Z - standardisation
    Z  = (Matrix{Float64}(X) .- self.call.means) ./self.call.std

    return predict_ind_sup(Z,ind_names,self.call.var_weights,self.svd.V)
end;

"""
    supvarPCA(self;...)

Supplementary variables (quantitative and/or qualitative) in Principal Components Analysis (PCA).

...
# Description
Performs the factor coordinates, square cosinus and square distance to origin of supplementary variables (quantitative and/or qualitative) with Principal Components Analysis (PCA).

# Required arguments
- `self`: a PCA NamedTuple

# Optional (default) arguments
- `X_quanti_sup`: a DataFrame of supplementary quantitative variables (default = nothing)
- `X_quali_sup`: a DataFrame of supplementary qualitative variables (default = nothing)
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is removed from DataFrame;
    - If `false`: the first columns is kepted in the DataFrame.

# Returns
A NamedTuple of NamedTuple containing all the results for the new variables including:
- `quanti_sup`: NamedTuple containing the results for the supplementary quantitative variables including:
    - `coord`: factor coordinates (scores) for the supplementary quantitative variables;
    - `cos2`: square cosinus for the supplementary quantitative variables;
    - `summary`: statistics (minimum, average, standard deviation, etc.) for the supplementary quantitative variables.

- `quali_sup`: NamedTuple containing the results for the supplementary qualitative variables including:
    - `coord`: factor coordinates (scores) for the supplementary categories
    - `cos2`: square cosinus for the supplementary categories
    - `vtest`: value-test for the supplementary categories
    - `dist2`: square distance to origin for the supplementary categories
    - `eta2`: square correlation ratio for the supplementary qualitative variables;
    - `summary`: count and proportions for the supplementary categories.

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Predict supplementary variables
julia> quanti_sup = get_dataset("decathlon",choice="quanti_sup");
julia> quali_sup = get_dataset("decathlon",choice="quali_sup"); 
julia> sup_var = supvarPCA(res_pca,X_quanti_sup=quanti_sup,X_quali_sup=quali_sup);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function supvarPCA(self::NamedTuple;
                   X_quanti_sup=nothing,
                   X_quali_sup=nothing,
                   first_col_as_index::Bool=true)::NamedTuple

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;

    # Check if PCA model
    if self.model != "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;   

    # Check if all are empty
    if X_quanti_sup === nothing && X_quali_sup === nothing  throw(ArgumentError("At least one shouldn't be empty."))  end;

    ##########################################################################################################
    ## Statistics for supplementary quantitatives variables
    ##########################################################################################################
    if X_quanti_sup !== nothing
        # Check if X_quanti_sup is a DataFrame
        if !isa(X_quanti_sup,AbstractDataFrame) throw(ArgumentError("X_quanti_sup must be a DataFrame.")) end;
        
        # Check if first columns is index names
        if first_col_as_index X_quanti_sup = X_quanti_sup[!,2:end] end;

        # Store all informations
        quanti_sup_ = predict_quanti_sup(X_quanti_sup,self.call.ind_weights,self.svd.U)
    else
        quanti_sup_ = nothing;
    end;

    ########################################################################################################
    ## Statistics for supplementary qualitatives variables
    ########################################################################################################
    if X_quali_sup !== nothing
        # Check if X_quali_sup is a DataFrame
        if !(X_quali_sup isa AbstractDataFrame) throw(ArgumentError("X_quali_sup must be a DataFrame.")) end;
        
        # Check if first columns is index names
        if first_col_as_index X_quali_sup = X_quali_sup[!,2:end] end;

        # Recode variables if two variables have at least one categories in common
        X_quali_sup = recode_cat_variable(X_quali_sup)
        
        # Statistics for qualitatives variables
        summary_quali_sup = freq_table(X_quali_sup)

        # Set supplementary categories names
        mod_sup_names = select(summary_quali_sup,:Modalities)

        # Conditional weighted average with original dataset
        barycentre = conditional_weighted_average(self.call.X[!,2:end],X_quali_sup,w=self.call.ind_weights)

        # Standardization
        Z_quali_sup = (Matrix{Float64}(barycentre[!,2:end]) .- self.call.means) ./ self.call.std

        # Supplementary categories factor coordinates
        quali_sup_coord = (Z_quali_sup .* transpose(self.call.var_weights)) * self.svd.V[:,1:self.call.n_components]

        # Supplementary categories squared distance to origin
        quali_sup_sqdisto = sum((Z_quali_sup .^2) .* transpose(self.call.var_weights),dims=2)
        
        # Supplementary categories cos2
        quali_sup_cos2 = (quali_sup_coord .^2) ./quali_sup_sqdisto

        # Supplementary categories values-tests
        n_rows, n_k = nrow(self.call.X), summary_quali_sup.Effectifs
        quali_sup_vtest = (quali_sup_coord ./transpose(sqrt.(self.eig[1:self.call.n_components,2]))) .* sqrt.(((n_rows - 1).*n_k)./(n_rows .- n_k))

        # square correlation ratio
        quali_sup_eta2 = reduce(vcat,[function_eta2(X_quali_sup,lab,self.ind.coord[!,2:end];w=self.call.ind_weights) for lab in names(X_quali_sup)])

        # Store all informations
        quali_sup_ = (; :barycentre => barycentre,
                        :coord => hcat(mod_sup_names,DataFrame(quali_sup_coord,self.call.dim_index)),
                        :cos2 => hcat(mod_sup_names,DataFrame(quali_sup_cos2,self.call.dim_index)),
                        :vtest => hcat(mod_sup_names,DataFrame(quali_sup_vtest,self.call.dim_index)),
                        :eta2 => quali_sup_eta2,
                        :dist2 => hcat(mod_sup_names,DataFrame("Sq. Dist" => vec(quali_sup_sqdisto))),
                        :summary => summary_quali_sup);
    else
        quali_sup_ = nothing;
    end;

    return (; :quanti_sup => quanti_sup_, :quali_sup => quali_sup_)
end;
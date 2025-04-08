"""
    FAMD(X;...)

Factor Analysis of Mixed Data (FAMD)

...
# Description
Performs Factor Analysis of Mixed Data (FAMD) with supplementary individuals, supplementary quantitative variables and supplementary categorical variables.

FAMD is a principal component method dedicated to explore data with both continuous and categorical variables. 
It can be seen roughly as a mixed between PCA and MCA. More precisely, 
the continuous variables are scaled to unit variance and the categorical variables are transformed 
into a disjunctive data table (crisp coding) and then scaled using the specific scaling of MCA. 
This ensures to balance the influence of both continous and categorical variables in the analysis. 
It means that both variables are on a equal foot to determine the dimensions of variability. 
This method allows one to study the similarities between individuals taking into account mixed 
variables and to study the relationships between all the variables.

# Details
FAMD includes standard Principal Component Analysis (PCA) and Multiple Correspondence Analysis (MCA) as special cases. If all variables are quantitative, standard PCA is performed. if all variables are qualitative, then standard MCA is performed.

# Required arguments
- `X`: DataFrame of shape (`n_samples`, `n_columns`)
    Training data, where `n_samples` is the numbeer of samples and `n_columns` is the number of columns.

# Optional (default) arguments
- `n_components`: an integer specified the number of dimensions kept in the results (by default 5).
- `ind_weights`: an optional individuals weights (by default, a vector of 1/(number of active individuals) for uniform individuals weights), the weights are given only for active individuals.
- `quanti_weights`: an optional quantitative variables weights (by default, a vector of 1 for uniform variables weights), the weights are given only for quantitative active variables.
- `quali_weights`: an optional qualitative variables weights (by default, a vector of 1 for uniform variables weights), the weights are given only for qualitative active variables.
- `ind_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary individuals.
- `quanti_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary quantitative variables.
- `quali_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary qualitative variables.
- `fill_missing`: a symbol specified method for missing treatment. Allowed values are:
    - `:mean` (default) for fill missing with mean
    - `:median` for fill missing with median
    - `drop` for dropping missing values.
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is used as index (individuals) names
    - If `false`: index (individuals) names are created.

A NamedTuple with following keys:
- `model`: string specifying the model fitted = 'famd'.
- `call`: namedtuple with some informations.
- `eig`: dataFrame containing all the eigenvalues, the difference between each eigenvalue, the percentage of variance and the cumulative percentage of variance.
- `svd`: namedtuple with results for generalied singular value decomposition (GSVD)
- `ind`: namedtuple of DataFrame containing all the results for the active individuals (factor coordinates, square cosinus, relative contributions, additionals informations).
- `quanti_var`: namedtuple of DataFrame containing all the results for the active quantitative variables (factor coordinates, square cosinus, relative contributions, additionals informations).
- `quali_var`: namedtuple of DataFrame containing all the results for the active qualitative variables (factor coordinates, square cosinus, relative contributions, additionals informations).
- `ind_sup`: namedtuple of DataFrame containing all the results for the supplementary individuals (factor coordinates, square cosinus,square distance to origin).
- `quanti_sup`: namedtuple of DataFrame containing all the results for the supplementative quantitative variables (factor coordinates, correlation between variables and axes, square cosinus).
- `quali_sup`: namedtuple of DataFrame containing all the results for supplementary qualitative variables (factor coordinates of each categories of each variables, vtest which is a criterion with a Normal distribution, and eta2 which is the square correlation coefficient between qualitative variable and an axe).
- `summary_quanti`: dataframe specifying summary statistics for quantitative variables (actives and/or supplementary).
- `summary_quali`: dataframe specifying summary statistics for qualitative variables (actives and/or supplementary).
- `association`: dataframe specifying association test. If qualitative are greater than 2.
- `chi2_test`: dataframe specifying chi2 contingency test. If qualitative are greater than 2.

# Examples
```julia-repl
julia> using Scientisttools
julia> autos2005 = get_dataset("autos2005");
julia> res_famd = FAMD(autos2005,ind_sup=39:45,quanti_sup=14:16,quali_sup=17);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

# References
- Escofier B, Pagès J (2023), Analyses Factorielles Simples et Multiples. 5ed, Dunod
- Husson F., Le S. and Pagès J. (2010). Exploratory Multivariate Analysis by Example Using R, Chapman and Hall.
- Husson F., Josse L, Lê S. & Mazet J. (2009). FactoMineR : Factor Analysis and Data Mining iwith R. R package version 2.11
- Lebart L., Piron M. & Morineau A. (2006). Statistique exploratoire multidimensionelle. Dunod Paris 4ed
- Lê, S., Josse, J., & Husson, F. (2008). FactoMineR: An R Package for Multivariate Analysis. Journal of Statistical Software, 25(1), 1–18. https://doi.org/10.18637/jss.v025.i01
- Pagès J. (2004). Analyse factorielle de donnees mixtes. Revue Statistique Appliquee. LII (4). pp. 93-111.
- Pagès J. (2013). Analyse factorielle multiple avec R : Pratique R. edp sciences
- Rakotomalala, Ricco (2020), Pratique des méthodes factorielles avec Python. Université Lumière Lyon 2, Version 1.0

# See Also
`get_famd_ind`, `get_famd_var`, `get_famd`, `dimdesc`, `predictFAMD`, `supvarFAMD`, `fviz_famd_ind`, `fviz_famd_col`, `fviz_famd_mod`, `fviz_famd_var`

...
"""
function FAMD(X::AbstractDataFrame;
              n_components::Int=5,
              ind_weights=nothing,
              quanti_weights=nothing,
              quali_weights=nothing,
              ind_sup=nothing,
              quanti_sup=nothing,
              quali_sup=nothing,
              fill_missing::Symbol=:mean,
              first_col_as_index::Bool=true)::NamedTuple
    
    ##################################################################################
    ## Check if X is a DataFrame
    ##################################################################################
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    ###################################################################################
    ## Check if supplementary individuals
    ####################################################################################
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
        throw(ArgumentError("Invalid fill missing name $(fill_missing)."))
    end;

    ############################################################################################
    ## Make a copy of the data
    ##############################################################################################
    Xtot = copy(X)

    ##############################################################################################
    ## Drop supplementary qualitatives variables
    ###############################################################################################
    if quali_sup !== nothing X = X[!,Not(quali_sup_label)] end;

    ###############################################################################################
    ## Drop supplementary quantitatives variables
    ###############################################################################################
    if quanti_sup !== nothing X = X[!,Not(quanti_sup_label)] end;
    
    ################################################################################################
    # Drop supplementary individuals
    ################################################################################################
    if ind_sup !== nothing
        # Extract supplementary individuals
        X_ind_sup = X[ind_sup_idx,:]
        X = X[Not(ind_sup_idx), :]
    end;

    ###############################################################################################
    # Set active individuals names
    ################################################################################################
    if first_col_as_index
        ind_names = DataFrame(Individuals=string.(X[!,1]))
        # Drop first columns
        X = X[!,2:end]
    else
        ind_names = DataFrame(Individuals=["ind" * string(i) for i in 1:nrow(X)])
    end;

    #
    # Factor Analysis of Mixed Data
    rec = recodevarfamd(X)

    # Extract all elements
    X, X_quanti, X_quali, dummies = rec.X, rec.quanti, rec.quali, rec.dummies
    n_rows, n_cont, n_cat, nb_moda = rec.n, rec.k1, rec.k2, rec.nb_moda

    ########################################################################################################
    ## Set individuals weights
    ########################################################################################################
    if ind_weights === nothing
        ind_weights = ones(n_rows)/n_rows
    elseif !isa(ind_weights,Array{<:Number,1})
        throw(ArgumentError("ind_weights must be a vector of number."))
    elseif length(ind_weights) !== n_rows
        throw(DimensionMismatch("ind_weights must be a vector with length $n_rows."))
    else
        ind_weights = ind_weights/sum(ind_weights)
    end;

    ##################################################################################################
    ## Set variables weights
    ##################################################################################################
    Z,means,stds,var_weights = ind_names, Vector{Float64}(),Vector{Float64}(),Vector{Float64}()

    #################################################################################################
    ## Set quantitative variables weights
    ##################################################################################################
    if n_cont > 0
        # Convert data to matrix
        vquant = Matrix{Float64}(X_quanti)
        # Statistics for quantitative variables
        summary_quanti = summary_stats(X_quanti)
        # Average and standard deviation
        wmean1, wstd1 = descrstatsw(vquant,weights(ind_weights));

        # Standardisation : Z = (X - mu)/sigma
        Z1 = DataFrame((vquant .- wmean1) ./ wstd1, names(X_quanti))

        # Concatenate
        Z, means, stds = hcat(Z,Z1), vcat(means,vec(wmean1)), vcat(stds,vec(wstd1))

        # Set quanti_weights
        if quanti_weights === nothing
            quanti_weights = ones(n_cont)
        elseif !isa(quanti_weights,Array{<:Number,1})
            throw(ArgumentError("quanti_weights must be a vector of number."))
        elseif length(quanti_weights) !== n_cont
            throw(DimensionMismatch("quanti_weights must be a vector with length $n_cont."))
        else
            quanti_weights = vec(quanti_weights)
        end;

        # Concatenate
        var_weights = vcat(var_weights,quanti_weights)
    end;

    #######################################################################################################
    ## Set qualitative variables weights
    #######################################################################################################
    if n_cat > 0
        # Summary statistics
        summary_quali = freq_table(X_quali)

        # Convert dummies to matrix
        dummies_mat = Matrix{Int}(dummies)
        # Normalize Z
        prop = sum(dummies_mat .* ind_weights,dims=1)
        # Compute weighted mean and weighted standard deviation
        wmean2, _ = descrstatsw(dummies_mat,weights(ind_weights))
        wstd2 = sqrt.(prop)

        # Creation of matrix Z
        Z2 = DataFrame((dummies_mat .- wmean2) ./ wstd2, names(dummies))
        
        # Concatenate
        Z, means, stds = hcat(Z,Z2), vcat(means,vec(wmean2)), vcat(stds,vec(wstd2))

        # Set qualitative variables weights
        if quali_weights === nothing
            quali_weights = ones(n_cat)
        elseif !isa(quali_weights,Array{<:Number,1})
            throw(ArgumentError("quali_weights must be a vector of number."))
        elseif length(quali_weights) != n_cat
            throw(DimensionMismatch("quali_weights must be a vector with length $n_cat."))
        else
            quali_weights = vec(quali_weights)
        end;

        # Extract nrow and proprow
        p_k = summary_quali.Proportions
        # Set modalities weights
        mod_weights = p_k .* reduce(vcat,[repeat([quali_weights[j]],nb_moda[j]) for j in 1:n_cat]);

        # Concatenate
        var_weights = vcat(var_weights,vec(mod_weights))
    end;

    # QR decomposition (to set maximum number of components)
    qr_dec = qr(Matrix(Z[!,2:end]))
    max_components = min(rank(Matrix(qr_dec.Q)),rank(Matrix(qr_dec.R)))
    
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

    # Store call informations

     # store call informations
     call_ = (; :Xtot => Xtot, 
                :X => hcat(ind_names,X), 
                :Z => Z,
                :ind_weights => ind_weights, 
                :var_weights => var_weights,
                :means => means, 
                :std => stds, 
                :n_components => n_components,
                :dim_index => dim_index,
                :ind_sup => ind_sup_idx,
                :quanti_sup => quanti_sup_label,
                :quali_sup => quali_sup_label,
                :rec => rec,
                :first_col_as_index => first_col_as_index)

    # Initialize named tuple
    res = (; :model => "famd", :call => call_);

    #################################################################################################
    ## Fit Principal Components Analysis (PCA)
    #################################################################################################
    fit_ = PCA(Z,standardize=false,n_components=max_components,ind_weights=ind_weights)

    #################################################################################################
    ## Statistics for active qualitative variables
    #################################################################################################
    if n_cat > 0
        # Concatenate with active qualitative variablles
        Z_quali = hcat(Z,X_quali)
        # Find index of qualitatives
        idx = [i for i in 1:ncol(Z_quali) if names(Z_quali)[i] in names(X_quali)]
        # Update fit
        fit_ = PCA(Z_quali,standardize=false,n_components=max_components,ind_weights=ind_weights,quali_sup=idx)

        # Extract elements
        quali_var_eta2 = fit_.quali_sup.eta2[!,1:(n_components+1)]

        # Store all informations
        quali_var_ = (; :coord => fit_.quali_sup.coord[!,1:(n_components+1)],
                        :contrib => fit_.var.contrib[(n_cont+1):end,1:(n_components+1)],
                        :cos2 => fit_.quali_sup.cos2[!,1:(n_components+1)],
                        :vtest => fit_.quali_sup.vtest[!,1:(n_components+1)],
                        :dist2 => fit_.quali_sup.dist2,
                        :infos => fit_.var.infos[(n_cont+1):end,:])
        
        # Others statistics
        if n_cat > 1
            association, chi2_test = fit_.association, fit_.chi2_test
        end;
    end;

    #########################################################################################
    ## Statistics for supplementary qualitatives variables
    #########################################################################################
    if quali_sup !== nothing
        # Select supplementary qualitatives variables
        X_quali_sup = Xtot[!,quali_sup_idx]

        # Transform to DataFrame if Vector
        if isa(X_quali_sup,AbstractVector) X_quali_sup = DataFrame(quali_sup_label => X_quali_sup) end;

        # Remove supplementary individuals
        if ind_sup !== nothing X_quali_sup = X_quali_sup[Not(ind_sup_idx),:] end;

        # Concatenate
        Z_quali_sup = leftjoin(Z,hcat(ind_names,X_quali_sup),on=:Individuals)
        # Find index of supplementative qualitative variables
        idx = [i for i in 1:ncol(Z_quali_sup) if names(Z_quali_sup)[i] in names(X_quali_sup)]
        # Update fit
        fit_ = PCA(Z_quali_sup,standardize=false,n_components=max_components,ind_weights=ind_weights,quali_sup=idx)

        # Store all informations
        quali_sup_ = (; :coord => fit_.quali_sup.coord[!,1:(n_components+1)],
                        :cos2 => fit_.quali_sup.cos2[!,1:(n_components+1)],
                        :vtest => fit_.quali_sup.vtest[!,1:(n_components+1)],
                        :eta2 => fit_.quali_sup.eta2[!,1:(n_components+1)],
                        :dist2 => fit_.quali_sup.dist2)

        # Store summary_quali_sup
        summary_quali_sup = fit_.summary_quali_sup

        # Update summary
        if n_cat > 0
            # Insert group columns
            insertcols!(summary_quali,1,:group => "active")
            insertcols!(summary_quali_sup,1,:group => "sup")

            # Concatenate
            summary_quali = vcat(summary_quali,summary_quali_sup)
        else
            res = @insert res.summary_quali_sup = summary_quali_sup;
        end;

        # Insert group
        if n_cat > 1
            insertcols!(chi2_test,1,:group => "active")
            insertcols!(association,1,:group => "active")
        end;

        # Update
        if ncol(X_quali_sup)>1
            association2, chi2_test2 = fit_.association, fit_.chi2_test
            insertcols!(chi2_test2,1,:group => "sup")
            insertcols!(association2,1,:group => "sup")

            if n_cat > 1
                # Concatenate DataFrame
                chi2_test,association = vcat(chi2_test,chi2_test2),vcat(association,association2)
            elseif n_cat === 1
                chi2_test, association = chi2_test2, association2
            else
                # Update NamedTuple
                res = @insert res.association = association2;
                res = @insert res.chi2_test = chi2_test2;
            end;
        end;

        # Others statistics
        if n_cat > 0
            ## Compute association test
            # Convert vector of tuple to vector of vector
            cols2 = [[col[1],col[2]] for col in collect(Iterators.product(names(X_quali_sup),names(X_quali)))]
            col_names2 = reduce(vcat,[DataFrame(Dict("variable1" => col[1],"variable2" => col[2])) for col in cols2])
            cont_table2 = Dict(col => freqtable(hcat(X_quali_sup,X_quali)[!,col],Symbol(col[1]),Symbol(col[2])) for col in cols2)
            # Chi2 statistic test
            chi2_test3 = hcat(col_names2,reduce(vcat,[DataFrame([chi2_contingency(cont_table2[col])[[:statistic,:dof,:pvalue]]]) for col in cols2]))
            association3 = hcat(col_names2,reduce(vcat,[association_measure(cont_table2[col])[!,2:end] for col in cols2]))

            insertcols!(chi2_test3,1,:group => "sup")
            insertcols!(association3,1,:group => "sup")

            # Insert group
            if n_cat > 1
                # Concatenate DataFrame
                chi2_test,association = vcat(chi2_test,chi2_test3),vcat(association,association3)
            else
                chi2_test, association = chi2_test3, association3
            end; 
        end;
        
        # Update namedtuple
        res = @insert res.quali_sup = quali_sup_;
    end;

    #########################################################################################
    ## Update Generalized Singular Value Decomposition (GSVD)
    #########################################################################################
    svd_ = (; :vs => fit_.svd.vs[1:max_components], :U => fit_.svd.U[:,1:n_components], :V => fit_.svd.V[:,1:n_components])
    
    # Update res NamedTuple
    res = @insert res.svd = svd_;

    #########################################################################################
    ## Eigenvalues
    #########################################################################################
    # Update res NamedTuple
    res = @insert res.eig = fit_.eig;

    #########################################################################################
    ## Statistics for active individuals
    #########################################################################################
    ind_ = (; :coord => fit_.ind.coord[!,1:(n_components+1)],
              :contrib => fit_.ind.contrib[!,1:(n_components+1)],
              :cos2 => fit_.ind.cos2[!,1:(n_components+1)],
              :infos => fit_.ind.infos)

    # Update NamedTuple
    res = @insert res.ind = ind_;

    #######################################################################################
    ## Statistics for active quantitative variables
    #######################################################################################
    if n_cont > 0
        # Contributions
        quanti_var_contrib = fit_.var.contrib[1:n_cont,1:(n_components+1)]
        # Square cosinus
        quanti_var_cos2 = fit_.var.cos2[1:n_cont,1:(n_components+1)]

        # Store all informations
        quanti_var_ = (; :coord => fit_.var.coord[1:n_cont,1:(n_components+1)], 
                         :contrib => quanti_var_contrib,
                         :cos2 => quanti_var_cos2, 
                         :infos => fit_.var.infos[1:n_cont,:])
        
        # Update res NamedTuple
        res = @insert res.quanti_var = quanti_var_;
    end;

    ####################################################################################
    ## Statistics for variables 
    ####################################################################################
    if n_cat > 0
        # Add to qualitative namedtuple
        if n_cont === 0
            # Update quali_var_ NamedTuple
            quali_var_ = @insert quali_var_.eta2 = quali_var_eta2;
        end;
        # Update res NamedTuple
        res = @insert res.quali_var = quali_var_;

        # Add variables informations if both qualitative and quantitative variables
        if n_cont > 0
            # Set variables names
            quanti_var_names, quali_var_names = select(quanti_var_contrib,:Variables), select(quali_var_eta2,:Variables)

            # Contribution des variables qualitatives
            quali_var_contrib = 100*(Matrix{Float64}(quali_var_eta2[!,2:end]) ./transpose(vec(fit_.eig[1:n_components,2])))
            # Cosinus carré des variables qualitatives
            quali_var_cos2 = (Matrix{Float64}(quali_var_eta2[!,2:end]) .^2) ./ vec(nb_moda .- 1)

            # Conversion in DataFrame
            quali_var_contrib = hcat(quali_var_names,DataFrame(quali_var_contrib,names(quali_var_eta2)[2:end]))
            quali_var_cos2 = hcat(quali_var_names,DataFrame(quali_var_cos2,names(quali_var_eta2)[2:end]))

            # Store all informations
            var_ = (; :coord => vcat(quanti_var_cos2,quali_var_eta2), 
                      :contrib => vcat(quanti_var_contrib,quali_var_contrib), 
                      :cos2 => vcat(hcat(quanti_var_names,DataFrame(Matrix{Float64}(quanti_var_cos2[!,2:end]).^2,names(quanti_var_cos2)[2:end])),quali_var_cos2))

            # Update res NamedTuple
            res = @insert res.var = var_;
        end;
    end;

    ###################################################################################################
    ## Statistics for supplementary individuals
    ###################################################################################################
    if ind_sup !== nothing
        # Number of supplementary individuals 
        n_rows_sup = nrow(X_ind_sup)
        # Set supplementary individuals names
        if first_col_as_index
            ind_sup_names = DataFrame(Individuals=string.(X_ind_sup[!,1]))
            # Drop first columns
            X_ind_sup = X_ind_sup[!,2:end] 
        else
            ind_sup_names = DataFrame(Individuals=["ind" * string(i) for i in (n_rows+1):(n_rows+n_rows_sup)])
        end;
        
        # Recode variables
        rec2 = recodevarfamd(X_ind_sup)

        # Extract elements
        X_ind_sup_quanti, X_ind_sup_quali, n_cont2, n_cat2 = rec2.quanti, rec2.quali, rec.k1, rec.k2

        # Initialize
        Z_ind_sup = ind_sup_names

        if n_cont2 > 0
            if n_cont !== n_cont2 throw(DimensionMismatch("The number of quantitative columns must be the same.")) end;
            # Standardise the data
            Z1_ind_sup = DataFrame((Matrix{Float64}(X_ind_sup_quanti) .- transpose(means[1:n_cont])) ./ transpose(stds[1:n_cont]),names(X_ind_sup_quanti))

            # Concatenate DataFrame
            Z_ind_sup = hcat(Z_ind_sup,Z1_ind_sup)
        end;

        if n_cat2 > 0
            if n_cat !== n_cat2 throw(DimensionMismatch("The number of qualitative columns must be the same.")) end;
            # Creation of supplementary dummies
            dummies_sup = zeros((n_rows_sup,ncol(dummies)));
            for i in 1:n_rows_sup
                values = [X_ind_sup_quali[i,j] for j in 1:n_cat]
                for k in 1:ncol(dummies)
                    if names(dummies)[k] in values
                        dummies_sup[i,k] = 1
                    end;
                end;
            end;

            # Standardiz data
            Z2_ind_sup = DataFrame((dummies_sup .- transpose(means[(n_cont+1):end]))./ transpose(stds[(n_cont+1):end]),names(dummies))

            # Concatenate DataFrame
            Z_ind_sup = hcat(Z_ind_sup,Z2_ind_sup)
        end;

        # Store all informations
        ind_sup_ = predict_ind_sup(Matrix{Float64}(Z_ind_sup[!,2:end]),ind_sup_names,ones(ncol(Z)-1),svd_.V)
 
        # Update namedtuple
        res = @insert res.ind_sup = ind_sup_;
    end;

    ##################################################################################
    ## Statistics for supplementary quantitatives variables
    ##################################################################################
    if quanti_sup !== nothing
        # Select supplementary quantitatives variables
        X_quanti_sup = Xtot[!,quanti_sup_idx]

        # Transform to DataFrame if Vector
        if isa(X_quanti_sup,AbstractVector) X_quanti_sup = DataFrame(quanti_sup_label => X_quanti_sup) end;

        # Remove supplementary individuals
        if ind_sup !== nothing X_quanti_sup = X_quanti_sup[Not(ind_sup_idx),:] end;

        # Store all informations
        quanti_sup_ = predict_quanti_sup(X_quanti_sup,ind_weights,svd_.U)
        
        # Extract elements
        summary_quanti_sup = quanti_sup_.summary

        # Insert group columns
        if n_cont === 0
             # Update namedtuple
            res = @insert res.summary_quanti_sup = summary_quanti_sup;
        elseif n_cont > 0
            insertcols!(summary_quanti,1, :group => "active")
            insertcols!(summary_quanti_sup,1, :group => "sup")
            # Concatenate
            summary_quanti = vcat(summary_quanti,summary_quanti_sup)
        end;
        
        # Update namedtuple
        res = @insert res.quanti_sup = quanti_sup_.quanti_sup;
    end;

    # Add summmary quanti
    if n_cont > 0 res = @insert res.summary_quanti = summary_quanti; end;

    # Add summary_quali
    if n_cat > 0
        # Update NamedTuple
        res = @insert res.summary_quali = summary_quali;
        if (n_cat > 1) || (quali_sup !== nothing)
            # Update NamedTuple
            res = @insert res.association = association;
            res = @insert res.chi2_test = chi2_test;
        end;
    end;

    return res
end;
"""
    predictFAMD(self,X;...)

Predict projection for new individuals with Factor Analysis of Mixed Data (FAMD).

...
# Description
Performs the coordinates, squared cosinus and square distance to origin of new individuals with Factor Analysis of Mixed Data (FAMD)

# Required arguments
- `self`: a FAMD NamedTuple
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
julia> using Scientisttools;
julia> autos2005 = get_dataset("autos2005");
julia> res_famd = FAMD(autos2005,ind_sup=39:45,quanti_sup=14:16,quali_sup=17);
julia> # Predict supplementary individuals
julia> donnee = get_dataset("autos2005",choice="ind_sup");
julia> ind_sup = predictFAMD(res_famd,afmd_ind_sup);
```

# Authors
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function predictFAMD(self::NamedTuple,
                     X::AbstractDataFrame;
                     first_col_as_index::Bool=true)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;

    # Check if X is a DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Check if FAMD model
    if self.model != "famd" throw(ArgumentError("self must be a FAMD NamedTuple.")) end;

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

    # Recode variables
    rec2 = recodevarfamd(X)

    # Extract elements
    X_quanti, X_quali, n_cont2, n_cat2 = rec2.quanti, rec2.quali, rec2.k1, rec2.k2

    # Active recoding
    n_cont, n_cat, dummies = self.call.rec.k1, self.call.rec.k2, self.call.rec.dummies

    # Initialize Z
    Z = ind_names

    if n_cont2 > 0
        if n_cont !== n_cont2 throw(DimensionMismatch("The number of quantitative columns must be the same.")) end;
        # Standardise the data
        Z1 = DataFrame((Matrix{Float64}(X_quanti) .- transpose(self.call.means[1:n_cont])) ./ transpose(self.call.std[1:n_cont]),names(X_quanti))
        # Concatenate DataFrame
        Z = hcat(Z,Z1)
    end;

    if n_cat2 > 0
        if n_cat !== n_cat2 throw(DimensionMismatch("The number of qualitative columns must be the same.")) end;
        # Creation of supplementary dummies
        dummies_sup = zeros((nrow(X),ncol(dummies)));
        for i in 1:nrow(X)
            values = [X_quali[i,j] for j in 1:n_cat]
            for k in 1:ncol(dummies)
                if names(dummies)[k] in values
                    dummies_sup[i,k] = 1
                end;
            end;
        end;
        # Standardize data
        Z2 = DataFrame((dummies_sup .- transpose(self.call.means[(n_cont+1):end]))./ transpose(self.call.std[(n_cont+1):end]),names(dummies))
        # Concatenate DataFrame
        Z = hcat(Z,Z2)
    end;

    # Store all informations
    return predict_ind_sup(Matrix{Float64}(Z[!,2:end]),ind_names,ones(ncol(Z)-1),self.svd.V)
end;
"""
    FAMD(self;...)

Supplementary variables in Factor Analysis of Mixed Data (FAMD)
...
# Description
Performs the coordinates, squared cosinus and squared distance to origin of supplementary variables with Factor Analysis of Mixed Data (FAMD).

# Required arguments
- `self`: a FAMD NamedTuple

# Optional arguments
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
julia> autos2005 = get_dataset("autos2005");
julia> res_famd = FAMD(autos2005,ind_sup=39:45,quanti_sup=14:16,quali_sup=17);
julia> # # Predict supplementary variables
julia> quanti_sup = get_dataset("autos2005",choice="quanti_sup");
julia> quali_sup = get_dataset("autos2005",choice="quali_sup");
julia> sup_var = supvarFAMD(res_famd,X_quanti_sup=quanti_sup,X_quali_sup=quali_sup);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function supvarFAMD(self::NamedTuple;
                    X_quanti_sup=nothing,
                    X_quali_sup=nothing,
                    first_col_as_index::Bool=true)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;

    # Check if FAMD model
    if self.model != "famd" throw(ArgumentError("self must be a FAMD NamedTuple.")) end;   

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
        if !isa(X_quali_sup,AbstractDataFrame) throw(ArgumentError("X_quali_sup must be a DataFrame.")) end;
        
        # Check if first columns is index names
        if first_col_as_index X_quali_sup = X_quali_sup[!,2:end] end;

        # Recode variables if two variables have at least one categories in common
        X_quali_sup = recode_cat_variable(X_quali_sup)
        
        # Statistics for qualitatives variables
        summary_quali_sup = freq_table(X_quali_sup)

        # Set supplementary categories names
        mod_sup_names = select(summary_quali_sup,:Modalities)

        # Set variables weights
        var_weights = ones(ncol(self.call.Z)-1)
        means, stds = mean(Matrix{Float64}(self.call.Z[!,2:end])), ones(ncol(self.call.Z)-1)

        # Conditional weighted average with original dataset
        barycentre = conditional_weighted_average(self.call.Z[!,2:end],X_quali_sup,w=self.call.ind_weights)

        # Standardization
        Z_quali_sup = (Matrix{Float64}(barycentre[!,2:end]) .- transpose(means)) ./ transpose(stds)

        # Supplementary categories factor coordinates
        quali_sup_coord = (Z_quali_sup .* transpose(var_weights)) * self.svd.V[:,1:self.call.n_components]

        # Supplementary categories squared distance to origin
        quali_sup_sqdisto = sum((Z_quali_sup .^2) .* transpose(var_weights),dims=2)
        
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
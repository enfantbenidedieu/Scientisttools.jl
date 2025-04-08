"""
    MCA(X;...)

Multiple Correspondence Analysis (MCA)

...
# Description

Performs Multiple Correspondence Analysis (MCA) over the data given in DataFrame X with supplementary individuals, supplementary categorical variables and supplementary quantitative variables. 

# Required arguments
- `X`: DataFrame of shape (n_samples, n_columns)
    Training data, where `n_samples` is the number of samples and `n_columns` is the number of columns.

# Optional (default) arguments
- `n_components`: an integer specified the number of dimensions kept in the results (by default 5).
- `ind_weights`: an optional individuals weights (by default, a vector of 1/(number of active individuals) for uniform individuals weights), the weights are given only for active individuals.
- `var_weights`: an optional variables weights (by default, a vector of 1/(number of active variables) for uniform variables weights), the weights are given only for active variables.
- `benzecri`: a boolean, default = true. If `true`, eigenvalues Benzecri correction is applied.
- `greenacre`: a boolean, default = true. If `true`, eigenvalues Greenacre correction is applied.
- `ind_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary individuals.
- `quali_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary qualitative variables.
- `quanti_sup`: an integer or a vector or an unitrange indicating the indexes of the supplementary quantitative variables.
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is used as index (individuals) names
    - If `false`: index (individuals) names are created.

# Returns
A NamedTuple with following keys:
- `model`: string specifying the model fitted = 'mca'.
- `call`: namedtuple with some informations.
- `eig`: dataFrame containing all the eigenvalues, the difference between each eigenvalue, the percentage of variance and the cumulative percentage of variance.
- `svd`: namedtuple with results for generalied singular value decomposition (GSVD)
- `ind`: namedtuple of DataFrame containing all the results for the active individuals (factor coordinates, square cosinus, relative contributions, additionals informations).
- `var`: namedtuple of DataFrame containing all the results for the active variables (factor coordinates, square cosinus, relative contributions, additionals informations).
- `others`: others informations for Multiple Correspondence Analysis (kaiser threshold, etc.).
- `ind_sup`: namedtuple of DataFrame containing all the results for the supplementary individuals (factor coordinates, square cosinus,square distance to origin).
- `quali_sup`: namedtuple of DataFrame containing all the results for supplementary qualitative variables (factor coordinates of each categories of each variables, vtest which is a criterion with a Normal distribution, and eta2 which is the square correlation coefficient between qualitative variable and an axe).
- `quanti_sup`: namedtuple of DataFrame containing all the results for the supplementative quantitative variables (factor coordinates, correlation between variables and axes, square cosinus).
- `summary_quanti_sup`: dataframe specifying summary statistics for supplementary quantitative variables.
- `summary_quali`: dataframe specifying summary statistics for qualitative variables (active and/or supplementary).
- `association`: dataframe specifying association test.
- `chi2_test`: dataframe specifying chi2 contingency test.

# Examples
```julia-repl
julia> using Scientisttools
julia> races_canines = get_dataset("races_canines");
julia> res_mca = MCA(races_canines,ind_sup=28:33,quali_sup=8,quanti_sup=9);
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
`get_mca_ind`, `get_mca_var`, `get_mca`, `dimdesc`, `reconst`, `predictMCA`, `supvarMCA`, `fviz_mca_ind`, `fviz_mca_var`, `fviz_mca_biplot`.

...
"""
function MCA(X::AbstractDataFrame;
             n_components::Int=5,
             ind_weights=nothing,
             var_weights=nothing,
             benzecri::Bool=true,
             greenacre::Bool=true,
             ind_sup=nothing,
             quali_sup=nothing,
             quanti_sup=nothing,
             first_col_as_index::Bool=true)::NamedTuple
    
    ##################################################################################
    ## Check if X is a DataFrame
    ##################################################################################
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    ##################################################################################
    ## Check if supplementary individuals
    ##################################################################################
    if ind_sup !== nothing
        if typeof(ind_sup) === Int
            ind_sup_idx = [Int(ind_sup)]
        elseif isa(ind_sup,AbstractVector)
            ind_sup_idx = reduce(vcat, ind_sup)
        elseif isa(ind_sup,AbstractUnitRange)
            ind_sup_idx = reduce(vcat,collect.([ind_sup]))
        end;
    else
        ind_sup_idx = nothing
    end;

    ##################################################################################
    ## Check if supplementary qualitatives variables
    #################################################################################
    if quali_sup !== nothing
        if typeof(quali_sup) === Int
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

    ###################################################################################
    ## Check if supplementary quantitatives variables
    ###################################################################################
    if quanti_sup !== nothing
        if typeof(quanti_sup) === Int
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

    ####################################################################################
    ## Make a copy of the data
    ####################################################################################
    Xtot = copy(X)

    ####################################################################################
    ## Drop supplementary quantitatives variables
    ####################################################################################
    if quanti_sup !== nothing
        X = X[!,Not(quanti_sup_label)]
    end;

    ####################################################################################
    ## Drop supplementary qualitatives variables
    ####################################################################################
    if quali_sup !== nothing
        X = X[!,Not(quali_sup_label)]
    end;

    ####################################################################################
    ## Drop supplementary individuals
    ####################################################################################
    if ind_sup !== nothing
        # Extract supplementary individuals
        X_ind_sup = X[ind_sup_idx,:]
        X = X[Not(ind_sup_idx), :]
    end;

    ####################################################################################
    # Create individuls and variables labels
    ####################################################################################
    # Create label for active individuals
    if first_col_as_index
        ind_names = DataFrame(Individuals=string.(X[!,1]))
        # Drop first columns
        X = X[!,2:end]
    else
        ind_names = DataFrame(Individuals=["ind" * string(i) for i in 1:nrow(X)])
    end;

    # Create label for active columns
    var_names = DataFrame(Variables=string.(names(X))) 

    # size of active elements
    n_rows, n_vars = nrow(ind_names), nrow(var_names)

    # Check if two variables have at least one categorie in common
    X = recode_cat_variable(X)

    # summary 
    summary_quali = freq_table(X);

    # Convert to matrix
    vqual = Matrix(X)

    ##########################################################################################################
    ## Compute association test
    ##########################################################################################################
    cols = collect(combinations(names(X),2))
    col_names = reduce(vcat,[DataFrame(Dict("variable1" => col[1],"variable2" => col[2])) for col in cols])
    cont_table = Dict(col => freqtable(X[!,col],Symbol(col[1]),Symbol(col[2])) for col in cols)

    # Chi2 statistic test
    chi2_test = hcat(col_names,reduce(vcat,[DataFrame([chi2_contingency(cont_table[col])[[:statistic,:dof,:pvalue]]]) for col in cols]))
    # Association test
    association = hcat(col_names,reduce(vcat,[association_measure(cont_table[col])[!,2:end] for col in cols]))

    # Dummies tables
    dummies = get_dummies(X);

    # Create label for actives modalities
    mod_names = DataFrame(Modalities=string.(names(dummies)))

    # Number of modalities
    n_mod = nrow(mod_names)

    # Convert dummies to matrix
    dummies_mat = Matrix{Int}(dummies)

    # Set individuals weights
    if ind_weights === nothing
        ind_weights = ones(n_rows)/n_rows
    elseif !(ind_weights isa Array{<:Number,1})
        throw(ArgumentError("ind_weights must be a vector of number."))
    elseif length(ind_weights) != n_rows
        throw(DimensionMismatch("ind_weights must be a vector with length $n_rows."))
    else
        ind_weights = ind_weights/sum(ind_weights)
    end;

     # set variables weights
     if var_weights === nothing
        var_weights = ones(n_vars)/n_vars
    elseif !(var_weights isa Array{<:Number,1})
        throw(ArgumentError("var_weights must be a vector of number."))
    elseif length(var_weights) != n_vars
        throw(DimensionMismatch("var_weights must be a vector with length $n_vars."))
    else
        var_weights = vec(var_weights)/sum(vec(var_weights))
    end;

    # Extract nrow and proprow
    I_k, p_k = summary_quali.Effectifs, summary_quali.Proportions

    # Creation of matrix Z (useful to SVD)
    Z = (dummies_mat ./ transpose(p_k)) .- 1

    # QR decomposition (to set maximum number of components)
    qr_dec = qr(Matrix(Z))
    max_components = Int(min(min(rank(Matrix(qr_dec.Q)),rank(Matrix(qr_dec.R))), n_rows - 1, n_mod - n_vars))

    ################################################################################################
    ## Set number of components
    ################################################################################################
    if n_components === nothing
        n_components = Int(max_components)
    elseif typeof(n_components) != Int
        throw(ArgumentError("n_components must be an integer."))
    elseif n_components < 1
        throw(ArgumentError("n_components must be equal or greater than 1."))
    else
        n_components = Int(min(n_components,max_components))
    end;

    # set columns
    dim_index = ["Dim." * string(i) for i in 1:n_components]

    ################################################################################################
    ## Set modalities weights
    ################################################################################################
    # Number of modalitie by qualitative variables
    nb_mod_var = [length(unique(X[!,j])) for j in 1:n_vars];
    # modalities weights
    mod_weights = p_k .* reduce(vcat,[repeat([var_weights[j]],nb_mod_var[j]) for j in 1:n_vars]);

    # store call informations
    call_ = (; :Xtot => Xtot, 
               :X => hcat(ind_names,X), 
               :dummies => hcat(ind_names,dummies),
               :Z => hcat(ind_names,DataFrame(Z,mod_names[!,1])),
               :ind_weights => ind_weights, 
               :var_weights => var_weights,
               :mod_weights => mod_weights,
               :n_components => n_components,
               :dim_index => dim_index,
               :ind_sup => ind_sup_idx,
               :quanti_sup => quanti_sup_label,
               :quali_sup => quali_sup_label,
               :first_col_as_index => first_col_as_index)

    # Initialize named tuple
    res = (; :model => "mca", :call => call_);

    ###############################################################################################
    ## Fit Factor Analysis model
    ###############################################################################################
    fit_ = fitFA(Z,max_components,n_components=n_components,row_weights=ind_weights,row_names=ind_names,col_weights=mod_weights,col_names=mod_names);

    # Extract elements and update namedtuple
    res = @insert res.svd = fit_.svd;
    res = @insert res.eig = fit_.eig;
    res = @insert res.ind = fit_.row;

    ## Others informations for variables
    # Normalized columns coordinates : see (Saporta, p235) or (Husson, p138)
    corrected_var_coord =  Matrix(fit_.col.coord[!,2:end]) .* transpose(fit_.svd.vs[1:n_components])

    # Variables value-test
    var_vtest = Matrix(fit_.col.coord[!,2:end]) .* sqrt.(((n_rows - 1).*I_k)./(n_rows .- I_k));
    
    # qualitative variables contributions
    quali_var_contrib = zeros((n_vars,n_components));
    for j in 1:n_vars
        contrib = Matrix{Float64}(filter(x -> x.Modalities in unique(vqual[:,j]),fit_.col.contrib)[:,2:end])
        quali_var_contrib[j,:] = sum(contrib,dims=1)
    end;

    # qualitative variables inertia
    quali_var_inertia =  (nb_mod_var .- 1)./n_vars
    quali_var_inertia_var_pct = 100*quali_var_inertia/sum(quali_var_inertia)

    # conversion en dataframe
    infos_quali_var = hcat(var_names,DataFrame([var_weights quali_var_inertia quali_var_inertia_var_pct],["Weight","Inertia","% Inertia"]));

    # qualitative variables square correlation ratio
    quali_var_eta2 = reduce(vcat,[function_eta2(X,lab,fit_.row.coord[!,2:end];w=ind_weights) for lab in names(X)])

    # Store all informations
    others_var_ = (;:corrected_coord => hcat(mod_names,DataFrame(corrected_var_coord,dim_index)),
                    :vtest => hcat(mod_names,DataFrame(var_vtest,dim_index)),
                    :var_contrib => hcat(var_names,DataFrame(quali_var_contrib,dim_index)),
                    :eta2 => quali_var_eta2,
                    :infos_var => infos_quali_var)

    # Update NamedTuple
    res = @insert res.var = merge(fit_.col,others_var_)

    # Eigenvalue threshold
    inertia = sum(quali_var_inertia)
    # Kaiser threshold
    kaiser_threshold = 1/n_vars
    kaiser_proportion_threshold = 100/inertia

    # Store others informations
    others_ = (; kaiser= (; threshold=kaiser_threshold, proportion_threshold=kaiser_proportion_threshold))

    # Update namedtuple
    res = @insert res.others = others_;

    ## Eigenvalue correction
    eigen_values = vec(fit_.eig[!,2])
    lambda = eigen_values[eigen_values .> kaiser_threshold]

    ################################################################################################
    ## Benzécri correction
    ################################################################################################
    if benzecri
        if length(lambda) > 0
            lambda_tilde = ((n_vars/(n_vars-1))*(lambda .- kaiser_threshold)).^2;
            prop_lambda_tilde = 100*lambda_tilde/sum(lambda_tilde);
            prop_lambda_tilde_cum = cumsum(prop_lambda_tilde);
            #convert to DataFrame
            benzecri_correction = DataFrame((; :Dimension => dim_index[1:length(lambda)],
                                               :Eigenvalue => lambda_tilde,
                                               :Proportion => prop_lambda_tilde,
                                               Symbol("Cum. proportion") => prop_lambda_tilde_cum))

            # Update namedtuple
            res = @insert res.benzecri_correction_ = benzecri_correction; 
        end;
    end;

    ###############################################################################################
    ## Greenacre correction
    ###############################################################################################
    if greenacre
        if length(lambda) > 0
            lambda_tilde = ((n_vars/(n_vars-1))*(lambda .- kaiser_threshold)).^2;
            s_tilde_tilde = (n_vars/(n_vars-1))*(sum(eigen_values .^2) - (n_mod - n_vars)/(n_vars .^2))
            tau = 100*(lambda_tilde/s_tilde_tilde);
            tau_cum = cumsum(tau);

            #convert to DataFrame
            greenacre_correction = DataFrame((; :Dimension => dim_index[1:length(lambda)],
                                                :Eigenvalue => lambda_tilde,
                                                :Proportion => tau,
                                                Symbol("Cum. proportion") => tau_cum))

            # Update namedtuple
            res = @insert res.greenacre_correction_ = greenacre_correction; 
        end;
    end;

    ###########################################################################################
    ## Statistics for supplementary individuals
    ###########################################################################################
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

        # Convert to matrix
        X_ind_sup_mat = Matrix{String}(X_ind_sup);
        # Creation of supplementary dummies
        dummies_sup = zeros((n_rows_sup,n_mod));
        for i in 1:n_rows_sup
            values = [X_ind_sup_mat[i,j] for j in 1:n_vars]
            for k in 1:n_mod
                if names(dummies)[k] in values
                    dummies_sup[i,k] = 1
                end;
            end;
        end;

        # Creation of matrix Z
        Z_ind_sup = (dummies_sup ./ transpose(p_k)) .- 1 

        # Strore all informations
        ind_sup_ = predict_ind_sup(Z_ind_sup,ind_sup_names,mod_weights,fit_.svd.V);

        # Update namedtuple
        res = @insert res.ind_sup = ind_sup_;
    end;

    #########################################################################
    ## Statistics for supplementary qualitative variables
    #########################################################################
    if quali_sup !== nothing
        # Select supplementary qualitatives variables
        X_quali_sup = Xtot[!,quali_sup_idx]

        # Transform to DataFrame if Vector
        if isa(X_quali_sup,AbstractVector) X_quali_sup = DataFrame(quali_sup_label => X_quali_sup) end;

        # Remove supplementary individuals
        if ind_sup !== nothing X_quali_sup = X_quali_sup[Not(ind_sup_idx),:] end;

        # Compute dummies tables
        dummies_sup = Matrix{Int}(get_dummies(X_quali_sup))

        # Extract nrow and proprow
        quali_sup_n_k, quali_sup_p_k = vec(sum(dummies_sup,dims=1)),vec(mean(dummies_sup,dims=1))

        # Creation of matrix Z
        Z_quali_sup = (dummies_sup ./ transpose(quali_sup_p_k)) .- 1

        # Supplementary qualitative categories
        quali_sup_coord = conditional_weighted_average(DataFrame(fit_.row.coord[!,2:end]),X_quali_sup,w=ind_weights)

        # Create supplementary modalities labels
        mod_sup_names = select(quali_sup_coord,:Modalities)

        # Correction with squared of eigenvalues
        quali_sup_coord = Matrix{Float64}(quali_sup_coord[!,2:end]) ./transpose(fit_.svd.vs[1:n_components])
        
        # Supplementary categories square distance to origin
        quali_sup_sqdisto = vec(sum((Z_quali_sup .^2) .* ind_weights,dims=1))
        
        # Supplementary categories square cosinus
        quali_sup_cos2 = (quali_sup_coord .^2) ./ quali_sup_sqdisto

        # Supplementary categories values-tests
        quali_sup_vtest = quali_sup_coord .* sqrt.(((n_rows - 1).*quali_sup_n_k)./(n_rows .- quali_sup_n_k))

        # Supplementary qualitative variables square correlation ratio
        quali_sup_eta2 = reduce(vcat,[function_eta2(X_quali_sup,lab,fit_.row.coord[!,2:end];w=ind_weights) for lab in names(X_quali_sup)])

        ## Compute association test
        # Convert vector of tuple to vector of vector
        cols2 = [[col[1],col[2]] for col in collect(Iterators.product(names(X_quali_sup),names(X)))]
        col_names2 = reduce(vcat,[DataFrame(Dict("variable1" => col[1],"variable2" => col[2])) for col in cols2])
        cont_table2 = Dict(col => freqtable(hcat(X_quali_sup,X)[!,col],Symbol(col[1]),Symbol(col[2])) for col in cols2)
        # Chi2 statistic test
        chi2_test2 = hcat(col_names2,reduce(vcat,[DataFrame([chi2_contingency(cont_table2[col])[[:statistic,:dof,:pvalue]]]) for col in cols2]))
        association2 = hcat(col_names2,reduce(vcat,[association_measure(cont_table2[col])[!,2:end] for col in cols2]))

        # Insert group
        insertcols!(chi2_test,1,:group => "active")
        insertcols!(chi2_test2,1,:group => "sup")
        insertcols!(association,1,:group => "active")
        insertcols!(association2,1,:group => "sup")

        # Concatenate DataFrame
        chi2_test,association = vcat(chi2_test,chi2_test2),vcat(association,association2)

        # Compute association test for more than one supplementary qualitative variable
        if length(names(X_quali_sup))>1
            # Compute association test
            cols3 = collect(combinations(names(X_quali_sup),2))
            col_names3 = reduce(vcat,[DataFrame(Dict("variable1" => col[1],"variable2" => col[2])) for col in cols3])
            cont_table3 = Dict(col => freqtable(X_quali_sup[!,col],Symbol(col[1]),Symbol(col[2])) for col in cols3)

            # Chi2 statistic test
            chi2_test3 = hcat(col_names3,reduce(vcat,[DataFrame([chi2_contingency(cont_table3[col])[[:statistic,:dof,:pvalue]]]) for col in cols3]))
            association3 = hcat(col_names3,reduce(vcat,[association_measure(cont_table3[col])[!,2:end] for col in cols3]))

            # Insert group column
            insertcols!(chi2_test3,1,:group => "sup")
            insertcols!(association3,1,:group => "sup")

            # Append
            chi2_test,association = vcat(chi2_test,chi2_test3),vcat(association,association3)
        end;

        # Statistics for qualitatives variables
        summary_quali_sup = freq_table(X_quali_sup)

        # Insert group columns
        insertcols!(summary_quali,1,:group => "active")
        insertcols!(summary_quali_sup,1,:group => "sup")

        # Concatenate
        summary_quali = vcat(summary_quali,summary_quali_sup)

        # Store all informations
        quali_sup_ = (; :coord => hcat(mod_sup_names,DataFrame(quali_sup_coord,dim_index)),
                        :cos2 => hcat(mod_sup_names,DataFrame(quali_sup_cos2,dim_index)),
                        :vtest => hcat(mod_sup_names,DataFrame(quali_sup_vtest,dim_index)),
                        :eta2 => quali_sup_eta2,
                        :dist2 => hcat(mod_sup_names,DataFrame("Sq. Dist" => vec(quali_sup_sqdisto))))

        # Update namedtuple
        res = @insert res.quali_sup = quali_sup_;
    end;

    ######################################################################################
    ## Statistics for supplementary quantitative variables
    ######################################################################################
    if quanti_sup !== nothing
        # Select supplementary quantitatives variables
        X_quanti_sup = Xtot[!,quanti_sup_idx]

        # Transform to DataFrame if Vector
        if isa(X_quanti_sup,AbstractVector) X_quanti_sup = DataFrame(quanti_sup_label => X_quanti_sup) end;

        # Remove supplementary individuals
        if ind_sup !== nothing X_quanti_sup = X_quanti_sup[Not(ind_sup_idx),:] end;

        # Store all informations
        quanti_sup_ = predict_quanti_sup(X_quanti_sup,ind_weights,fit_.svd.U)

        # Update namedtuple
        res = @insert res.quanti_sup = quanti_sup_.quanti_sup;
        res = @insert res.summary_quanti_sup = quanti_sup_.summary;
    end;

    # Update NamedTuple
    res = @insert res.chi2_test = chi2_test;
    res = @insert res.association = association;
    res = @insert res.summary_quali = summary_quali;

    return res
end;

"""
    predictMCA(self,X;...)

Predict projection for new individuals with Multiple Correspondence Analysis (MCA).

...
# Description
Performs the coordinates, square cosinus and square distance to origin of new individuals with Multiple Correspondence Analysis (MCA).

# Required arguments
- `self`: a MCA NamedTuple
- `X`: a DataFrame in which to look for variables with which to predict. X must contain columns with the same names as the original data.

# Optionals arguments
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
julia> races_canines = get_dataset("races_canines");
julia> res_mca = MCA(races_canines,ind_sup=28:33,quali_sup=8,quanti_sup=9);
julia> # Predict supplementary individuals
julia> donnee = get_dataset("races_canines",choice="ind_sup");
julia> ind_sup = predictMCA(res_mca,donnee);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function predictMCA(self::NamedTuple,
                    X::AbstractDataFrame;
                    first_col_as_index::Bool=true)
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;

    # Check if X is a DataFrame
    if !isa(X,AbstractDataFrame) throw(ArgumentError("X must be a DataFrame.")) end;

    # Check if MCA model
    if self.model != "mca" throw(ArgumentError("self must be a MCA NamedTuple.")) end;

    # Set number
    n_rows, n_rows_sup, n_mod = nrow(self.call.X), nrow(X), nrow(self.var.coord)

    # Active modalities names
    mod_names = names(self.call.dummies)[2:end]

    # Set individuals names
    if first_col_as_index
        ind_names = DataFrame(Individus=string.(X[!,1]))
        # Drop first columns
        X = X[!,2:end] 
    else
        ind_names = DataFrame(Individus=["ind" * string(i) for i in (n_rows+1):(n_rows+n_rows_sup)])
    end;

    # Check if columns re aligned
    if ncol(X) != ncol(self.call.X[!,2:end]) throw(ArgumentError("Columns aren't aligned.")) end;

    # Creation of dummies data
    dummies = zeros((n_rows_sup,n_mod));
    for i in 1:n_rows_sup
        values = [Matrix{String}(X)[i,j] for j in 1:ncol(X)]
        for k in 1:n_mod
            if mod_names[k] in values
                dummies[i,k] = 1
            end;
        end;
    end;

    # Creation of matrix Z
    Z = (dummies ./ transpose(self.summary_quali.Proportions[1:n_mod])) .- 1 

    return predict_ind_sup(Z,ind_names,self.call.mod_weights,self.svd.V)
end;

"""
    supvarMCA(self;...)

Supplementary variables (quantitative and/or qualitative) in Multiple Correspondence Analysis (MCA).

...
# Description
Performs the factor coordinates, square cosinus and square distance to origin of supplementary variables (quantitative and/or qualitative) with Multiple Correspondence Analysis (MCA).

# Required arguments
- `self`: a MCA NamedTuple

# Optional (default) arguments
- `X_quali_sup`: a DataFrame of supplementary qualitative variables (default = nothing)
- `X_quanti_sup`: a DataFrame of supplementary quantitative variables (default = nothing)
- `first_col_as_index`: a boolean, default = true
    - If `true`: the first columns is removed from DataFrame;
    - If `false`: the first columns is keep in the DataFrame.

# Returns
A NamedTuple of NamedTuple containing all the results for the new variables including:
- `quali_sup`: NamedTuple containing the results for the supplementary qualitative variables including:
    - `coord`: factor coordinates (scores) for the supplementary categories;
    - `cos2`: square cosinus for the supplementary categories;
    - `vtest`: value-test for the supplementary categories;
    - `dist2`: square distance to origin for the supplementary categories;
    - `eta2`: square correlation ratio for the supplementary qualitative variables;
    - `summary`: count and proportions for the supplementary categories.

- `quanti_sup`: NamedTuple containing the results for the supplementary quantitative variables including:
    - `coord`: factor coordinates (scores) for the supplementary quantitative variables;
    - `cos2`: square cosinus for the supplementary quantitative variables;
    - `summary`: statistics (minimum, average, standard deviation, etc.) for the supplementary quantitative variables.

# Examples
```julia-repl
julia> using Scientisttools
julia> races_canines = get_dataset("races_canines");
julia> res_mca = MCA(races_canines,ind_sup=28:33,quali_sup=8,quanti_sup=9);
julia> # Predict supplementary variables
julia> quali_sup = get_dataset("races_canines",choice="quali_sup");
julia> quanti_sup = get_dataset("races_canines",choice="quanti_sup");
julia> sup_var = supvarMCA(res_mca,X_quali_sup=quali_sup,X_quanti_sup=quanti_sup);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function supvarMCA(self::NamedTuple;
                   X_quali_sup=nothing,
                   X_quanti_sup=nothing,
                   first_col_as_index::Bool=true)

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;

    # Check if MCA model
    if self.model != "mca" throw(ArgumentError("self must be a MCA NamedTuple.")) end;   

    # Check if all are empty
    if X_quanti_sup === nothing && X_quali_sup === nothing  throw(ArgumentError("At least one shouldn't be empty."))  end;
    
    ###########################################################################################
    ## Statistics for supplementary qualitative variables
    ###########################################################################################
    if X_quali_sup !== nothing
        # Check if X_quali_sup is a DataFrame
        if !isa(X_quali_sup,AbstractDataFrame) throw(ArgumentError("X_quali_sup must be a DataFrame.")) end;
        
        # Check if first columns is a vector of string
        if first_col_as_index X_quali_sup = X_quali_sup[!,2:end] end;
        
        # Recode variables if two variables have at least one categories in common
        X_quali_sup = recode_cat_variable(X_quali_sup)

        # Set number of active rows
        n_rows = nrow(self.call.X)

        # Conditional weighted average with factor coordinates
        quali_sup_coord = conditional_weighted_average(self.ind.coord[!,2:end],X_quali_sup,w=self.call.ind_weights)

        # Create modalities labels
        mod_sup_names = select(quali_sup_coord,:Modalities)

        # Correction with squared of eigenvalues
        quali_sup_coord = Matrix{Float64}(quali_sup_coord[!,2:end]) ./ transpose(sqrt.(self.eig[1:self.call.n_components,2]))

        # Compute dummies tables
        dummies = Matrix{Int}(get_dummies(X_quali_sup))

        # Extract nrow and proprow
        n_k, p_k = vec(sum(dummies,dims=1)), vec(mean(dummies,dims=1))

        # Creation of matrix Z
        Z_quali_sup = (dummies ./ transpose(p_k)) .- 1

        # Square distance to origin
        quali_sup_sqdisto = vec(sum((Z_quali_sup .^2) .* self.call.ind_weights,dims=1))
        
        # Square cosinus
        quali_sup_cos2 = (quali_sup_coord .^2) ./ quali_sup_sqdisto

        # Values-tests
        quali_sup_vtest = quali_sup_coord .* sqrt.(((n_rows - 1).*n_k)./(n_rows .- n_k))

        # Square correlation ratio
        quali_sup_eta2 = reduce(vcat,[function_eta2(X_quali_sup,lab,self.ind.coord[!,2:end];w=self.call.ind_weights) for lab in names(X_quali_sup)])

        # Statistics for qualitatives variables
        summary_quali_sup = freq_table(X_quali_sup)

        # Store all informations
        quali_sup_ = (; :coord => hcat(mod_sup_names,DataFrame(quali_sup_coord,self.call.dim_index)),
                        :cos2 => hcat(mod_sup_names,DataFrame(quali_sup_cos2,self.call.dim_index)),
                        :vtest => hcat(mod_sup_names,DataFrame(quali_sup_vtest,self.call.dim_index)),
                        :eta2 => quali_sup_eta2,
                        :dist2 => hcat(mod_sup_names,DataFrame("Sq. Dist" => quali_sup_sqdisto)),
                        :summary => summary_quali_sup);
    else
        quali_sup_ = nothing;
    end;

    #########################################################################################
    ## Statistics for supplementary quantitative variables
    #########################################################################################
    if X_quanti_sup !== nothing
        # Check if X_quanti_sup is a DataFrame
        if !isa(X_quanti_sup,AbstractDataFrame) throw(ArgumentError("X_quanti_sup must be a DataFrame.")) end;
        
        # Check if first columns is a vector of string
        if first_col_as_index X_quanti_sup = X_quanti_sup[!,2:end] end;
        
        # Store all informations
        quanti_sup_ = predict_quanti_sup(X_quanti_sup,self.call.ind_weights,self.svd.U)
    else
        quanti_sup_ = nothing;
    end;

    return (; :quali_sup => quali_sup_, :quanti_sup => quanti_sup_)
end;
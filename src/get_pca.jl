"""
    get_pca_ind(self)

Extract the results for individuals - PCA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active individuals from Principal Components Analysis (PCA) outputs.

# Arguments
- `self`: a PCA NamedTuple

# Returns
A NamedTuple of DataFrame containing all the results for the active individuals including:
- `coord`: factor coordinates (scores) of the individuals;
- `cos2`: square cosinus of the individuals;
- `contrib`: relative contributions of the individuals;
- `infos`: additionals informations (weight, square distance to origin, inertia and percentage of inertia) of the individuals.

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Extract the results for individuals
julia> ind = get_pca_ind(res_pca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_pca_ind(self::NamedTuple)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;
    return self.ind
end;

"""
    get_pca_var(self)

Extract the results for variables - PCA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active variables from Principal Components Analysis (PCA) outputs.

# Arguments
- `self`: a PCA NamedTuple

# Returns
A NamedTuple of DataFrame containing all the results for the active variables including:
- `coord`: factor coordinates (scores) of the variables;
- `cos2`: square cosinus of the variables;
- `contrib`: relative contributions of the variables;
- `infos`: additionals informations (weight, square distance to origin, inertia and percentage of inertia) of the variables.

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Extract the results for variables
julia> ind = get_pca_var(res_pca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_pca_var(self::NamedTuple)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;
    return self.var
end;

"""
    get_pca(self;...)

Extract the results for individuals/variables - PCA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active individuals/variables from Principal Components Analysis (PCA) outputs.
- `get_pca()`: Extract the results for individuals and variables
- `get_pca_ind()`: Extract the results for individuals only
- `get_pca_var()`: Extract the results for variables onfly

# Required arguments
- `self`: a PCA NamedTuple

# Optional arguments
- `choice`: the element to subset from the output. Allowed valued are:
    - `"ind"` for individuals
    - `"var"` for variables

# Returns
A NamedTuple of DataFrame containing all the results for the active individuals/variables including:
- `coord`: factor coordinates (scores) of the individuals/variables;
- `cos2`: square cosinus of the individuals/variables;
- `contrib`: relative contributions of the individuals/variables;
- `infos`: additionals informations (weight, square distance to origin, inertia and percentage of inertia) of the individuals/variables.

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Extract the results for individuals
julia> ind = get_pca(res_pca) # or get_pca(res_pca,choice="ind");
julia> # Extract the results for variables
julia> var = get_pca(res_pca,choice="var");
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_pca(self::NamedTuple;choice::AbstractString="ind")::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;
    # Check if choice in vector
    if !(choice in ["ind","var"]) throw(ArgumentError("'choice' should be one of 'ind' or 'var'.")) end;

    if choice == "ind"
        return get_pca_ind(self)
    else
        return get_pca_var(self)
    end;
end;
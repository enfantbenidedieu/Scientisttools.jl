"""
    get_mca_ind(self)

Extract the results for individuals - MCA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active individuals from Multiple Correspondence Analysis (PCA) outputs.

# Arguments
- `self`: a MCA NamedTuple

# Returns
A NamedTuple of DataFrame containing all the results for the active individuals including:
- `coord`: factor coordinates (scores) of the individuals;
- `cos2`: square cosinus of the individuals;
- `contrib`: relative contributions of the individuals;
- `infos`: additionals informations (weight, square distance to origin, inertia and percentage of inertia) of the individuals.

# Examples
```julia-repl
julia> using Scientisttools
julia> races_canines = get_dataset("races_canines");
julia> res_mca = MCA(races_canines,ind_sup=28:33,quali_sup=8,quanti_sup=9);
julia> # Extract individuals informations
julia> ind = get_mca_ind(res_mca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_mca_ind(self::NamedTuple)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if MCA model
    if self.model !== "mca" throw(ArgumentError("self must be a MCA NamedTuple.")) end;
    return self.ind
end;

"""
    get_mca_var(self)

Extract the results for variables - MCA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active variables from Multiple Correspondence Analysis (MCA) outputs.

# Arguments
- `self`: a MCA NamedTuple

# Returns
A NamedTuple of DataFrame containing all the results for the active variables including:
- `coord`: factor coordinates (scores) of the variables;
- `cos2`: square cosinus of the variables;
- `contrib`: relative contributions of the variables;
- `infos`: additionals informations (weight, square distance to origin, inertia and percentage of inertia) of the variables.

# Examples
```julia-repl
julia> using Scientisttools
julia> races_canines = get_dataset("races_canines");
julia> res_mca = MCA(races_canines,ind_sup=28:33,quali_sup=8,quanti_sup=9);
julia> # Extract variables informations
julia> vars = get_mca_var(res_mca);
```

# Authors
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_mca_var(self::NamedTuple)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if MCA model
    if self.model !== "mca" throw(ArgumentError("self must be a MCA NamedTuple.")) end;
    return self.var
end;

"""
    get_mca(self;...)

Extract the results for individuals/variables - MCA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active individuals/variables from Multiple Correspondence Analysis (MCA) outputs.
- `get_mca()`: Extract the results for individuals and variables
- `get_mca_ind()`: Extract the results for individuals only
- `get_mca_var()`: Extract the results for variables onfly

# Required arguments
- `self`: a MCA NamedTuple

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
julia> races_canines = get_dataset("races_canines","all");
julia> res_mca = MCA(races_canines,ind_sup=28:33,quali_sup=8,quanti_sup=9);
julia> # Extract individuals informations
julia> ind = get_mca(res_mca) # or get_mca(res_mca,choice="ind");
julia> # Extract variables informations
julia> vars = get_mca(res_mca,choice="var");
```

# Authors
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_mca(self::NamedTuple;choice::AbstractString="ind")::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if MCA model
    if self.model !== "mca" throw(ArgumentError("self must be a MCA NamedTuple.")) end;
    # Check if choice in vector
    if !(choice in ["ind","var"]) throw(ArgumentError("'choice' should be one of 'ind' or 'var'.")) end;

    if choice == "ind"
        return get_mca_ind(self)
    else
        return get_mca_var(self)
    end;
end;
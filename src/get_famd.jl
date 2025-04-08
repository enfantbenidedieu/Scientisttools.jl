"""
    get_famd_ind(self)

Extract the results for individuals - FAMD

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active individuals from Factor Analysis of Mixed Data (FAMD) outputs.

# Arguments
- `self`: a FAMD NamedTuple

# Returns
A NamedTuple of DataFrame containing all the results for the active individuals including:
- `coord`: factor coordinates (scores) of the individuals;
- `cos2`: square cosinus of the individuals;
- `contrib`: relative contributions of the individuals;
- `infos`: additionals informations (weight, square distance to origin, inertia and percentage of inertia) of the individuals.

# Examples
```julia-repl
julia> using Scientisttools
julia> autos2005 = get_dataset("autos2005");
julia> res_famd = FAMD(autos2005,ind_sup=39:45,quanti_sup=14:16,quali_sup=17);
julia> # Extract the results for individuals
julia> ind = get_famd_ind(res_famd);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_famd_ind(self::NamedTuple)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if FAMD model
    if self.model != "famd" throw(ArgumentError("self must be a FAMD NamedTuple.")) end;
    return self.ind
end;

"""
    get_famd_var(self;...)

Extract the results for variables - FAMD

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active variables from Factor Analysis of Mixed Data (FAMD) outputs.

# Required arguments
- `self`: a FAMD NamedTuple

" Optional arguments
- `choice`: the element to subset from the output. Allowed valued are:
    - `"var"` (default) for active variables
    - `"quanti_var"` for active quantitative variables
    - `"quali_var"` for active qualitative variables

# Returns
A NamedTuple of DataFrame containing all the results for the active variables including:
- `coord`: factor coordinates (scores) of the variables;
- `cos2`: square cosinus of the variables;
- `contrib`: relative contributions of the variables;

# Examples
```julia-repl
julia> using Scientisttools
julia> autos2005 = get_dataset("autos2005");
julia> res_famd = FAMD(autos2005,ind_sup=39:45,quanti_sup=14:16,quali_sup=17);
julia> # Extract the results for quantitative variables
julia> quanti_var = get_famd_var(res_famd,choice="quanti_var");
julia> # Extract the results for qualitative variables
julia> quali_var = get_famd_var(res_famd,choice="quanti_var");
julia> # Extract the results for variables
julia> vars = get_famd_var(res_famd) # or get_famd_var(res_famd,choice="var");
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_famd_var(self::NamedTuple;choice::AbstractString="var")::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if FAMD model
    if self.model != "famd" throw(ArgumentError("self must be a FAMD NamedTuple.")) end;
    # Check if choice in vector
    if !(choice in ["var","quanti_var","quali_var"]) throw(ArgumentError("'choice' should be one of 'var', 'quanti_var','quali_var'.")) end;

    if choice === "quali_var"
        if !(:quali_var in keys(self)) throw(ArgumentError("No qualitative variables.")) end;
        return self.quali_var
    end;

    if choice === "quanti_var"
        if !(:quanti_var in keys(self)) throw(ArgumentError("No quantitative variables.")) end;
        return self.quali_var
    end;

    if choice === "var"
        if !(:var in keys(self)) throw(ArgumentError("No mixed data.")) end;
        return self.var
    end;
end;

"""
    get_famd(self;...)

Extract the results for individuals/variables - FAMD

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active individuals/variables from Factor Analysis of Mixed Data (FAMD) outputs.
- `get_famd()`: Extract the results for individuals and variables
- `get_famd_ind()`: Extract the results for individuals only
- `get_famd_var()`: Extract the results for variables onfly

# Required arguments
- `self`: a FAMD NamedTuple

# Optional arguments
- `choice`: the element to subset from the output. Allowed valued are:
    - `"ind"` (default) for active individuals
    - `"var"` for active variables
    - `"quanti_var"` for active quantitative variables
    - `"quali_var"` for active qualitative variables

# Returns
A NamedTuple of DataFrame containing all the results for the active individuals/variables including:
- `coord`: factor coordinates (scores) of the individuals/variables;
- `cos2`: square cosinus of the individuals/variables;
- `contrib`: relative contributions of the individuals/variables;

# Examples
```julia-repl
julia> using Scientisttools
julia> autos2005 = get_dataset("autos2005");
julia> res_famd = FAMD(autos2005,ind_sup=39:45,quanti_sup=14:16,quali_sup=17);
julia> # Extract the results for individuals
julia> ind = get_famd(res_famd) # or get_famd(res_famd,choice="ind") ;
julia> # Extract the results for quantitative variables
julia> quanti_var = get_famd(res_famd,choice="quanti_var");
julia> # Extract the results for qualitative variables
julia> quali_var = get_famd(res_famd,choice="quanti_var");
julia> # Extract the results for variables
julia> vars = get_famd(res_famd,choice="var");

```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_famd(self::NamedTuple;choice::AbstractString="ind")::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if FAMD model
    if self.model != "famd" throw(ArgumentError("self must be a FAMD NamedTuple.")) end;
    # Check if choice in vector
    if !(choice in ["ind","var","quanti_var","quali_var"]) throw(ArgumentError("'choice' should be one of 'ind', 'var', 'quanti_var','quali_var'.")) end;

    if choice == "ind"
        return get_famd_ind(self)
    else 
        return get_famd_var(self,choice=choice)
    end;
end;
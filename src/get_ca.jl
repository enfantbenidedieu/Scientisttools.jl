"""
    get_ca_row(self)

Extract the results for rows - CA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active rows from Correspondence Analysis (CA) outputs.

# Arguments
- `self`: a CA NamedTuple

# Returns
A NamedTuple of DataFrame containing all the results for the active rows including:
- `coord`: factor coordinates (scores) of the rows;
- `cos2`: square cosinus of the rows;
- `contrib`: relative contributions of the rows;
- `infos`: additionals informations (weight, margin, square distance to origin, inertia and percentage of inertia) of the rows.

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Extract rows informations
julia> rows = get_ca_row(res_ca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_ca_row(self::NamedTuple)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if CA model
    if self.model !== "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;
    return self.row
end;

"""
    get_ca_col(self)

Extract the results for columns - CA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active columns from Correspondence Analysis (CA) outputs.

# Arguments
- `self`: a CA NamedTuple

# Returns
A NamedTuple of DataFrame containing all the results for the active columns including:
- `coord`: factor coordinates (scores) of the columns;
- `cos2`: square cosinus of the columns;
- `contrib`: relative contributions of the columns;
- `infos`: additionals informations (margin, square distance to origin, inertia and percentage of inertia) of the columns.

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Extract columns informations
julia> cols = get_ca_col(res_ca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_ca_col(self::NamedTuple)::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if CA model
    if self.model !== "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;
    return self.col
end;

"""
    get_ca(self;...)

Extract the results for rows/columns - CA

...
# Description
Extract all the results (factor coordinates, square cosinus, relative contributions and others informations) of the active rows/columns from Correspondence Analysis (CA) outputs.
- `get_ca()`: Extract the results for rows and columns
- `get_ca_row()`: Extract the results for rows only
- `get_ca_col()`: Extract the results for columns onfly

# Required arguments
- `self`: a CA NamedTuple

# Optional arguments
- `choice`: the element to subset from the output. Allowed valued are:
    - `"row"` for rows
    - `"col"` for columns

# Returns
A NamedTuple of DataFrame containing all the results for the active rows/columns including:
- `coord`: factor coordinates (scores) of the rows/variables;
- `cos2`: square cosinus of the rows/columns;
- `contrib`: relative contributions of the rows/columns;
- `infos`: additionals informations (margin, square distance to origin, inertia and percentage of inertia) of the rows/columns.

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Extract rows informations
julia> rows = get_ca(res_ca) # or get_ca(res_ca,choice="row");
julia> # Extract columns informations
julia> cols = get_ca(res_ca,choice="col");
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_ca(self::NamedTuple;choice::AbstractString="row")::NamedTuple
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if CA model
    if self.model !== "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;
    # Check if choice in vector
    if !(choice in ["row","col"]) throw(ArgumentError("'choice' should be one of 'row' or 'col'.")) end;
    
    if choice == "row"
        return get_ca_row(self)
    else
        return get_ca_col(self)
    end;
end;
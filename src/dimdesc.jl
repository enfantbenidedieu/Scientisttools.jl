"""
    dimdecs(self,axis;...)

Dimension description

...
# Description
This function is designed to point out the variables and the categories that are the most characteristic according to each dimension obtained by a Factor Analysis.

# Required arguments
- `self`: a PCA, CA, MCA and FAMD NamedTuple

# Optional (default) arguments
- `axis`: an integer or a vector or an unitrange indicating the indexes of axis
- `proba`: the significance threshold considered to characterized the category (by default 0.05)

# Returns
a NamedTuple containing
- `quanti`: the description of the dimensions by the quantitative variables. The variables are sorted.
- `quali`: the description of the dimensions by the qualitative variables. The variables are sorted.
- `category`: the description of dimensions by categories of qualitative variables.

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Dimension description
julia> dim_desc = dimdesc(res_pca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function dimdesc(self::NamedTuple;
                 axis=nothing,
                 proba::AbstractFloat=0.05)::NamedTuple

    function desc(x::AbstractDataFrame,y::AbstractVector,w,proba)
        # Split the data
        x_quanti = splitmix(x).quanti
        x_quali = splitmix(x).quali

        # initialize
        res = (;)
        if x_quanti !== nothing
            quanti = contdesc(x_quanti,y,w=w,proba=proba)
            if nrow(quanti) > 0 res = @insert res.quanti = quanti end;
        end;

        if x_quali !== nothing
            cat = catdesc(x_quali,y,proba=proba)
            if nrow(cat.quali) > 0 res = @insert res.quali = cat.quali end;
            if nrow(cat.category) > 0 res = @insert res.category = cat.category end;
        end;
        return res
    end;

    # Check if self is an object of NamedTuple
    if !(self isa NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if model in vector
    if !(self.model in ["pca","ca","mca","famd"]) throw(ArgumentError("self must be a PCA, CA, MCA or FAMD NamedTuple.")) end;

    if axis !== nothing
        if isa(axis,Int)
            axis = [Int(axis)]
        elseif isa(axis,AbstractVector)
            axis = reduce(vcat, axis)
        elseif isa(axis,AbstractUnitRange)
            axis = reduce(vcat,collect.([axis]))
        end;
    end;
    
    if self.model === "ca"
        # Extract rows and columns coordinates
        row_coord, col_coord = self.row.coord, self.col.coord

        # Add supplementary
        if :row_sup in keys(self) row_coord = vcat(row_coord, self.row_sup.coord) end;
        if :col_sup in keys(self) col_coord = vcat(col_coord, self.col_sup.coord) end;

        # Set row names and columns names
        row_names, col_names = row_coord[!, [1]], col_coord[!, [1]]
        # Update 
        row_coord, col_coord = row_coord[!,2:end], col_coord[!,2:end]

        if axis !== nothing
            row_coord, col_coord = row_coord[!,axis], col_coord[!,axis]
            if isa(row_coord,AbstractVector) row_coord = DataFrame(self.call.dim_index[axis] => row_coord) end;
            if isa(col_coord,AbstractVector) col_coord = DataFrame(self.call.dim_index[axis] => col_coord) end;
        end;

        function corrdim(idx)
            # Rows
            rows = row_coord[!,idx]
            if isa(rows,AbstractVector) rows = DataFrame(idx => rows) end;
            rename!(rows,Symbol(idx) => :coord)
            rows = sort(hcat(row_names,rows),:coord, rev=true)

            # Columns
            cols = col_coord[!,idx]
            if isa(cols,AbstractVector) cols = DataFrame(idx => cols) end;
            rename!(cols,Symbol(idx) => :coord)
            cols = sort(hcat(col_names,cols),:coord, rev=true)
            return (; :row => rows, :col => cols)
        end;
        corr_dim = NamedTuple(Symbol(idx) => corrdim(idx) for idx in names(row_coord));
    else
        # Extract global data
        X = self.call.Xtot
        # Remove first columns if rownames
        if self.call.first_col_as_index X = X[!,2:end] end;
        # Remove supplementary individuals
        if self.call.ind_sup !== nothing X = X[Not(self.call.ind_sup),2:end] end;

        # Extract individuals factor coordinates and weights
        ind_coord, ind_weights = self.ind.coord[!,2:end], self.call.ind_weights

        if axis !== nothing
            ind_coord = ind_coord[!,axis] 
            if isa(ind_coord,AbstractVector) ind_coord = DataFrame(self.call.dim_index[axis] => ind_coord) end;
        end;
        corr_dim = NamedTuple(Symbol(idx) => desc(X,ind_coord[!,idx],ind_weights,proba) for idx in names(ind_coord));
    end;

    return corr_dim
end;
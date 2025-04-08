"""
    get_eig(self)

Extract the eigenvalues/variances of dimensions

...
# Description
Eigenvalues correspond to the amount of the variation explained by each principal component.

# Arguments
- `self`: a PCA, CA, MCA and FAMD NamedTuple

# Returns
A DataFrame containing eigenvalue, difference, variance percent and cumulative variance of percent

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Extract eigenvalues
julia> eig = get_eig(res_pca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_eig(self::NamedTuple)::DataFrame
    # Check if self is an object of NamedTuple
    if !(self isa NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if model in vector
    if !(self.model in ["pca","ca","mca","famd"]) throw(ArgumentError("self must be a PCA, CA, MCA or FAMD NamedTuple.")) end;
    return self.eig
end;

"""
    get_eigenvalue(self)

Extract the eigenvalues/variances of dimensions

...
See `get_eig`

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function get_eigenvalue(self::NamedTuple)
    return get_eig(self)
end

"""
    fviz_screeplot(self;....)

Visualize the eigenvalues/proportions/cumulative of dimensions

...
# Description
This function suppport the results of multiple general factor analysis methods such as PCA (Principal Components Analysis), CA (Correspondence Analysis), MCA (Multiple Correspondence Analysis), etc...

# Required arguments
- `self`: a PCA, CA, MCA or FAMD NamedTuple

# Optional arguments
- `choice`: a symbol specified the data to plotted. Allowed values are `:proportion`, `:eigenvalue` or  `:cumulative`.
- `geom_type`: a symbol specifying the geometry to be used for the graph. Allowed values are `:bar` for bar plot, `:line` for lineplot or `[:bar, :line]` to use both types.
- `y_lim`: a numeric list of length 2 specifying the range of the plotted 'Y' values (by default = nothing).
- `bar_fill`: fill color for bar plot (by default = "steelblue").
- `bar_color`: outline color for bar plot (by default = "steelblue").
- `line_color`: color for line plot (by default = "black").
- `line_style`: a symbol specifying the linestyle of the line (by default = :dash)
- `ncp`: a numeric value specifying the number of dimensions to be shown.
- `add_labels`: a boolean (by default = false). If true, labels are added at the top of bars or points showing the information retained by each dimension.
- `label_color`: color for labels (by default = "black").
- `title`: a string corresponding to the title of the graph you draw (by default = nothing and a title is chosen).
- `x_label`: a string specifying the label text of x (by default = nothing and a x_label is chosen).
- `y_label`: a string specifying the label text of y (by default = nothing and a y_label is chosen).
- `hide_spine`: a boolean, default = false. If true, spines are hidden.
- `figsize`: the aspect of figure (by default = (800,500)).

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Extract eigenvalues
julia> eig = get_eig(res_pca);
julia> # Visualize eigenvalue
julia> p = fviz_screplot(res_pca);
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fviz_screeplot(self::NamedTuple;
                        choice::Symbol=:proportion,
                        geom=[:bar,:line],
                        y_lim = nothing,
                        bar_fill::Symbol=:steelblue,
                        bar_color::Symbol=:steelblue,
                        line_color::Symbol=:black,
                        line_style::Symbol=:dash,
                        ncp=10,
                        add_labels::Bool=false,
                        label_color::Symbol=:black,
                        title=nothing,
                        x_label=nothing,
                        y_label=nothing,
                        hide_spine::Bool=false,
                        figsize=(800,500))

    # Extract eigenvalue
    eig = get_eigenvalue(self)[1:min(ncp,self.call.n_components),:]

    # Organize dat
    if choice == :eigenvalue
        eig = vec(eig[!,"Eigenvalue"])
        text_labels = string.(round.(eig,digits=3))
        if y_label === nothing y_label = "Eigenvalue" end;
    elseif choice == :proportion
        eig = vec(eig[!,"Proportion"]) /100
        text_labels = string.(round.(100*eig,digits=2)) .* "%"
    elseif choice == :cumulative
        eig = vec(eig[!,"Cum. Proportion"]) /100
        text_labels = string.(round.(100*eig,digits=2)) .* "%"
        if y_label === nothing y_label = "Cumulative % of explained variances" end;
    else
        throw(ArgumentError("'choice' should be one of ':eigenvalue', ':proportion', ':cumulative'."))
    end;

    #
    if isa(geom,AbstractVector)
        intersect = [x for x in geom if x in [:bar,:line]]
        if length(intersect) == 0 throw(ArgumentError("Specified value(s) for the argument geom aren't allowed.")) end;
    elseif Base.isidentifier(geom)
        if !(geom in [:bar,:line]) throw(ArgumentError("Invalid geom name $(geom)")) end;
        geom = [geom]
    end;

    # Initialize
    fig = Figure(size=figsize)

    # Set title
    if title === nothing title = "Scree plot" end;
    # Set x labels
    if x_label === nothing x_label = "Dimensions" end;
    # Set y labels
    if y_label === nothing y_label = "% of explained variances" end;

    xticks_values,xticks_labels = collect(1:length(eig)),[string(x) for x in 1:length(eig)]
    ax = Axis(fig[1, 1],title = title,xlabel = x_label,ylabel = y_label,xticks = (xticks_values, xticks_labels))

    # Add bar
    if :bar in geom barplot!(ax,1:length(eig),eig,color = bar_color) end;

    # Add lines and scatter
    if :line in geom 
        lines!(ax,1:length(eig),eig,color = line_color,linestyle = line_style)
        scatter!(ax,1:length(eig),eig,color = line_color)
    end;

    # Add labels
    if add_labels text!(ax,1:length(eig),eig,text=text_labels,color=label_color,position=(:center,:top)) end;
    
    # Set y limits
    if y_lim !== nothing ylims!(ax,low = y_lim[1],high=y_lim[2]) end
    # Hide spine
    if hide_spine hidespines!(ax) end;

    return fig
end;
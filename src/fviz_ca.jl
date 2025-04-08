"""
    fviz_ca_row(self;...)

Visualize Correspondence Analysis (CA) - Graph of rows

...
# Description
Correspondence analysis (CA) is an extension of Principal Component Analysis (PCA) suited to analyze frequencies formed by two categorical variables. `fviz_ca_row` provides CairoMakie based elegant visualization of CA outputs from Julia functions.

# Required arguments
- `self`: a CA NamedTuple

# Optional arguments
- `axis`: a numeric vector/tuple of length 2 specifying the dimension to be plotted (by default =[0,1]).
- `x_lim`: a numeric vector/tuple of length 2 specifying the range of plotted 'x' values (by default = nothing).
- `y_lim`: a numeric vector/tuple of length 2 specifying the range of plotted 'y' values (by default = nothing).
- `x_label`: a string specifying the label text of x (by default = nothing).
- `y_label`: a string specifying the label text of y (by default = nothing).
- `title`: a string corresponding to the title of the graph you draw (by default = nothing).
- `color`: a symbol specifying the color for active individuals (by default = :black)
- `geom`: a symbol specifying the geometry to be used for the graph. Allowed values are the combinaison of [:point,:text]. Use ':point'  (to show only points); ':text' to show only labels; [:point,:text] to show both types.
- `point_size`: a numeric value specifying the marker size (by default = 10).
- `text_size`: a numeric value specifying the label size (by default = 15).
- `marker`: a symbol specifying the marker style (by default = :circle). For more, see https://docs.makie.org/stable/reference/plots/scatter#Default-markers
- `palette`: a symbol specifying the palette color of the colorbar (by default = :Hiroshige).
- `add_color_bar`: a boolean, default = true. If `true`, the colorbar is plotted, else the colorbar isn't plotted.
- `legend_title`: a string corresponding to the title of the legend (by default = nothing).
- `habillage`: a string corresponding the name of qualitative variable (by default = nothing).
- `row_sup`: a boolean to either add or not supplementary rows (by default = true).
- `row_sup_color`: a symbol specifying a color for the supplementary rows points (by default = :blue).
- `row_sup_marker`: a marker style for the supplementary rows points (by default = :rect).
- `row_sup_point_size`: a numeric value specifying the supplementary rows marker size (by default = 10).
- `row_sup_text_size`: a numeric value specifying the supplementary rows label size (by default = 15).
- `quali_sup`: a boolean to either add or not supplementary categories (by default = true).
- `quali_sup_color`: a symbol specifying a color for the supplementary categories points (by default = :red).
- `quali_sup_marker`: a marker style for the supplementary categories points (by default = :rect).
- `quali_sup_point_size`: a numeric value specifying the supplementary categories marker size (by default = 10).
- `quali_sup_text_size`: a numeric value specifying the supplementary categories label size (by default = 15).
- `add_hline`: a boolean to either add or not a horizontal ligne (by default = true).
- `hline_color`: a symbol specifying the horizontal ligne color (by default = :black).
- `hline_style`: a symbol specifying the horizontal ligne style (by default = :dash). Allowed values are : ':solid', ':dash', ':dashdot' or ':dashdotdot'.
- `hline_width`: a numeric specifying the horizontal line width (by default = 2).
- `add_vline`: a boolean to either add or not a vertical ligne (by default = true).
- `vline_color`: a symbol specifying the vertical ligne color (by default = :black).
- `vline_style`: a symbol specifying the vertical ligne style (by default = :dash). Allowed values are : ':solid', ':dash', ':dashdot' or ':dashdotdot'.
- `vline_width`: a numeric specifying the vertical line width (by default = 2).
- `ha`: horizontal alignment (by default = :left). Allowed values are : ':left', ':center' or ':right'.
- `va`: vertical alignment (by default = :center). Allowed values are : ':bottom', ':baseline', ':center' or ':top'.
- `hide_spine`: a boolean, default = false. If true, spines are hidden.
- `figsize`: the aspect of figure (by default = (800,500)).

# Return
a CairoMakie plot

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Graph of rows
julia> fig = fviz_ca_row(res_ca);
julia> fig
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fviz_ca_row(self::NamedTuple;
                      axis::AbstractVector=[1,2],
                      x_lim=nothing,
                      y_lim=nothing,
                      x_label=nothing,
                      y_label=nothing,
                      title = nothing,
                      color=:black,
                      geom=[:point,:text],
                      point_size=10,
                      text_size=15,
                      marker::Symbol=:circle,
                      palette::Symbol=:Hiroshige,
                      add_color_bar::Bool=true,
                      legend_title=nothing,
                      habillage = nothing,
                      row_sup::Bool=true,
                      row_sup_color::Symbol=:blue,
                      row_sup_marker::Symbol=:rect,
                      row_sup_point_size=10,
                      row_sup_text_size=15,
                      quali_sup::Bool=true,
                      quali_sup_color::Symbol=:red,
                      quali_sup_marker::Symbol=:hexagon,
                      quali_sup_point_size=10,
                      quali_sup_text_size=15,
                      add_hline::Bool=true,
                      hline_color::Symbol=:black,
                      hline_style::Symbol=:dash,
                      hline_width=2,
                      add_vline::Bool=true,
                      vline_color::Symbol=:black,
                      vline_style::Symbol=:dash,
                      vline_width=2,
                      ha::Symbol=:right,
                      va::Symbol=:bottom,
                      hide_spine::Bool=false,
                      figsize=(800,500))

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if CA model
    if self.model !== "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;

    if (length(axis) !== 2) || (axis[1] < 1) || (axis[2] > self.call.n_components) || (axis[1] > axis[2]) throw(ArgumentError("You must pass a valid 'axis'.")) end;

    # Check geom
    if isa(geom,AbstractVector)
        intersect = [x for x in geom if x in [:point,:text]]
        if length(intersect) == 0 throw(ArgumentError("Specified value(s) for the argument geom aren't allowed.")) end;
    elseif Base.isidentifier(geom)
        if !(geom in [:point,:text]) throw(ArgumentError("Invalid geom_type name $(geom)")) end;
        geom = [geom]
    end;

    # Add active data
    X = convert_to_float(self.call.X[!,2:end])
    coord = hcat(self.row.coord,X)

    # Add supplementary columns
    if :col_sup in keys(self)
        X_col_sup = self.call.Xtot[!,self.call.col_sup]
        if isa(X_col_sup,AbstractVector) X_col_sup = DataFrame(self.call.col_sup => X_col_sup) end;
        if self.call.row_sup !== nothing X_col_sup = X_col_sup[Not(self.call.row_sup),:] end;
        X_col_sup = convert_to_float(X_col_sup)
        coord = hcat(coord,X_col_sup)
    end;

    # Add supplementary quantitative variables
    if :quanti_sup in keys(self)
        X_quanti_sup = self.call.Xtot[!,self.call.quanti_sup]
        if isa(X_quanti_sup,AbstractVector) X_quanti_sup = DataFrame(self.call.quanti_sup => X_quanti_sup) end;
        if self.call.row_sup !== nothing X_quanti_sup = X_quanti_sup[Not(self.call.row_sup),:] end;
        X_quanti_sup = convert_to_float(X_quanti_sup)
        coord = hcat(coord,X_quanti_sup)
    end;

    # Add supplementary quantitative variables
    if :quali_sup in keys(self)
        X_quali_sup = self.call.Xtot[!,self.call.quali_sup]
        if isa(X_quali_sup,AbstractVector) X_quali_sup = DataFrame(self.call.quali_sup => X_quali_sup) end;
        if self.call.row_sup !== nothing X_quali_sup = X_quali_sup[Not(self.call.row_sup),:] end;
        coord = hcat(coord,X_quali_sup)
    end;

    if isa(color,String)
        if color === "cos2" 
            color_vec = vec(sum(Matrix{Float64}(self.row.cos2[!,2:end])[:,axis],dims=2)) 
            if legend_title === nothing legend_title = "Cos2" end;
        end;

        if color === "contrib" 
            color_vec = vec(sum(Matrix{Float64}(self.row.contrib[!,2:end])[:,axis],dims=2))
            if legend_title === nothing legend_title = "Contrib" end;
        end;

        if color in names(coord)[2:end]
            color_vec = vec(coord[!,color])
            if !isa(color_vec,Array{<:Number,1}) throw(ArgumentError("'color' must be a numeric variable.")) end;
            if legend_title === nothing legend_title = color end;
        end;
    elseif isa(color,Array{<:Number,1})
        color_vec = vec(color)
        if length(color_vec) !== nrow(coord) throw(ArgumentError("'color' must be a numeric vector of length $(nrow(coord)).")) end;
        if legend_title === nothing legend_title = "Cont_var" end;
    end;

    # For habillage
    if isa(habillage,String)
        if !(habillage in names(coord)[2:end]) throw(ArgumentError("$habillage not in DataFrame.")) end;
        vsqual = vec(coord[!,habillage])
        if legend_title === nothing legend_title = habillage end;
    elseif isa(habillage,Array{String,1})
        vsqual = vec(habillage)
        if length(vsqual) !== nrow(coord) throw(ArgumentError("'habillage' must be a string vector of length $(nrow(coord)).")) end;
        if legend_title === nothing legend_title = "Categorie" end;
    end;

    # Set x and y
    x_dim, y_dim = axis[1]+1, axis[2]+1

    # Set x labels
    if x_label === nothing x_label = "Dim." * string(axis[1]) *" ("*string(round.(self.eig[axis[1],4],digits=2))*"%)" end;
    # Set y labels
    if y_label === nothing y_label = "Dim." * string(axis[2]) *" ("*string(round.(self.eig[axis[2],4],digits=2))*"%)" end;
    # Set title
    if title === nothing title = "Rows Factor Map - CA" end;

    fig = Figure(size=figsize)
    ax = Axis(fig[1,1],xlabel=x_label,ylabel=y_label,title=title)

    if habillage === nothing
        if (isa(color,String) && (color in reduce(vcat,[["cos2","contrib"],names(coord)[2:end]]))) || isa(color,Array{<:Number,1})
            color_range = [minimum(color_vec),maximum(color_vec)]
            if :point in geom scatter!(ax,coord[!,x_dim],coord[!,y_dim],marker=marker,color=color_vec,colormap=palette,colorrange = color_range,markersize=point_size) end;
            if :text in geom text!(ax,coord[!,x_dim],coord[!,y_dim],text=coord[!,1],color=color_vec,colormap=palette,colorrange = color_range,fontsize=text_size,position=(ha,va)) end;
            if add_color_bar Colorbar(fig[1, 2],colorrange = color_range,colormap = palette, label = legend_title) end;
        else
            if :point in geom scatter!(ax,coord[!,x_dim],coord[!,y_dim],marker=marker,color=color,markersize=point_size) end;
            if :text in geom text!(ax,coord[!,x_dim],coord[!,y_dim],text=coord[!,1],color=color,fontsize=text_size,position=(ha,va)) end;
        end;
    else
        uq_group,colors = sort(unique(vsqual)),[:blue,:red,:green,:steelblue]
        Random.seed!(1234);
        color_group = Dict(zip(uq_group,sample(colors, length(uq_group); replace=false)))
        for g in uq_group
            idx = [i for i in 1:length(vsqual) if vsqual[i]==g]
            if :point in geom scatter!(ax,coord[idx,x_dim],coord[idx,y_dim],label=g,color=color_group[g],marker=marker,markersize=point_size) end;
            if :text in geom text!(ax,coord[idx,x_dim],coord[idx,y_dim],text=coord[idx,1],color=color_group[g],fontsize=text_size,position=(ha,va)) end;
        end;
        axislegend(ax,habillage)
    end;

    # Add supplementary rows coordinates
    if row_sup
        if :row_sup in keys(self) 
            row_sup_coord = self.row_sup.coord
            if :point in geom scatter!(ax,row_sup_coord[!,x_dim],row_sup_coord[!,y_dim],marker=row_sup_marker,color=row_sup_color,markersize=row_sup_point_size)end;
            if :text in geom text!(ax,row_sup_coord[!,x_dim],row_sup_coord[!,y_dim],text=row_sup_coord[!,1],color=row_sup_color,fontsize=row_sup_text_size,position=(ha,va)) end;
        end;
    end;

    # Add supplementary qualitatives variables
    if quali_sup
        if :quali_sup in keys(self)
            if habillage === nothing
                quali_sup_coord = self.quali_sup.coord
                if :point in geom scatter!(ax,quali_sup_coord[!,x_dim],quali_sup_coord[!,y_dim],marker=quali_sup_marker,color=quali_sup_color,markersize=quali_sup_point_size)end;
                if :text in geom text!(ax,quali_sup_coord[!,x_dim],quali_sup_coord[!,y_dim],text=quali_sup_coord[!,1],color=quali_sup_color,fontsize=quali_sup_text_size,position=(ha,va)) end;
            end;
        end;
    end;

    # Set x limits
    if x_lim !== nothing xlims!(ax,low = x_lim[1],high=x_lim[2]) end;
    # Set y limits
    if y_lim !== nothing ylims!(ax,low = y_lim[1],high=y_lim[2]) end;

    # Add horizontal line (hline)
    if add_hline hlines!(ax,0,color = hline_color,linestyle = hline_style,linewidth=hline_width) end;
    # Add vertical line (vline)
    if add_vline vlines!(ax,0,color = vline_color,linestyle = vline_style,linewidth=vline_width) end;
    # Hide spine
    if hide_spine hidespines!(ax) end;

    return fig
end;
"""
    fviz_ca_col(self;...)

Visualize Correspondence Analysis (CA) - Graph of columns

...
# Description
Correspondence analysis (CA) is an extension of Principal Component Analysis (PCA) suited to analyze frequencies formed by two categorical variables. `fviz_ca_col` provides CairoMakie based elegant visualization of CA outputs from Julia functions.

# Required arguments
- `self`: a CA NamedTuple

# Optional arguments
- `axis`: a numeric vector/tuple of length 2 specifying the dimension to be plotted (by default =[0,1]).
- `x_lim`: a numeric vector/tuple of length 2 specifying the range of plotted 'x' values (by default = nothing).
- `y_lim`: a numeric vector/tuple of length 2 specifying the range of plotted 'y' values (by default = nothing).
- `x_label`: a string specifying the label text of x (by default = nothing).
- `y_label`: a string specifying the label text of y (by default = nothing).
- `title`: a string corresponding to the title of the graph you draw (by default = nothing).
- `color`: a symbol specifying the color for active individuals (by default = :black)
- `geom`: a symbol specifying the geometry to be used for the graph. Allowed values are the combinaison of [:point,:text]. Use ':point'  (to show only points); ':text' to show only labels; [:point,:text] to show both types.
- `point_size`: a numeric value specifying the marker size (by default = 10).
- `text_size`: a numeric value specifying the label size (by default = 15).
- `marker`: a symbol specifying the marker style (by default = :circle). For more, see https://docs.makie.org/stable/reference/plots/scatter#Default-markers
- `palette`: a symbol specifying the palette color of the colorbar (by default = :Hiroshige).
- `add_color_bar`: a boolean, default = true. If `true`, the colorbar is plotted, else the colorbar isn't plotted.
- `legend_title`: a string corresponding to the title of the legend (by default = nothing).
- `col_sup`: a boolean to either add or not supplementary columns (by default = true).
- `col_sup_color`: a symbol specifying a color for the supplementary columns points (by default = :blue).
- `col_sup_marker`: a marker style for the supplementary columns points (by default = :rect).
- `col_sup_point_size`: a numeric value specifying the supplementary columns marker size (by default = 10).
- `col_sup_text_size`: a numeric value specifying the supplementary columns label size (by default = 15).
- `add_hline`: a boolean to either add or not a horizontal ligne (by default = true).
- `hline_color`: a symbol specifying the horizontal ligne color (by default = :black).
- `hline_style`: a symbol specifying the horizontal ligne style (by default = :dash). Allowed values are : ':solid', ':dash', ':dashdot' or ':dashdotdot'.
- `hline_width`: a numeric specifying the horizontal line width (by default = 2).
- `add_vline`: a boolean to either add or not a vertical ligne (by default = true).
- `vline_color`: a symbol specifying the vertical ligne color (by default = :black).
- `vline_style`: a symbol specifying the vertical ligne style (by default = :dash). Allowed values are : ':solid', ':dash', ':dashdot' or ':dashdotdot'.
- `vline_width`: a numeric specifying the vertical line width (by default = 2).
- `ha`: horizontal alignment (by default = :left). Allowed values are : ':left', ':center' or ':right'.
- `va`: vertical alignment (by default = :center). Allowed values are : ':bottom', ':baseline', ':center' or ':top'.
- `hide_spine`: a boolean, default = false. If true, spines are hidden.
- `figsize`: the aspect of figure (by default = (800,500)).

# Return
a CairoMakie plot

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Graph of columns
julia> fig = fviz_ca_col(res_ca);
julia> fig
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fviz_ca_col(self::NamedTuple;
                      axis::AbstractVector=[1,2],
                      x_lim=nothing,
                      y_lim=nothing,
                      x_label=nothing,
                      y_label=nothing,
                      title = nothing,
                      color=:black,
                      geom=[:point,:text],
                      point_size=10,
                      text_size=15,
                      marker::Symbol=:circle,
                      palette::Symbol=:Hiroshige,
                      add_color_bar::Bool=true,
                      legend_title=nothing,
                      col_sup::Bool=true,
                      col_sup_color::Symbol=:blue,
                      col_sup_marker::Symbol=:rect,
                      col_sup_point_size=10,
                      col_sup_text_size=15,
                      add_hline::Bool=true,
                      hline_color::Symbol=:black,
                      hline_style::Symbol=:dash,
                      hline_width=2,
                      add_vline::Bool=true,
                      vline_color::Symbol=:black,
                      vline_style::Symbol=:dash,
                      vline_width=2,
                      ha::Symbol=:right,
                      va::Symbol=:bottom,
                      hide_spine::Bool=false,
                      figsize=(800,500))

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if CA model
    if self.model !== "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;
    # Check axis
    if (length(axis) !== 2) || (axis[1] < 1) || (axis[2] > self.call.n_components) || (axis[1] > axis[2]) throw(ArgumentError("You must pass a valid 'axis'.")) end;

    # Check geom
    if isa(geom,AbstractVector)
        intersect = [x for x in geom if x in [:point,:text]]
        if length(intersect) == 0 throw(ArgumentError("Specified value(s) for the argument geom aren't allowed.")) end;
    elseif Base.isidentifier(geom)
        if !(geom in [:point,:text]) throw(ArgumentError("Invalid geom_type name $(geom)")) end;
        geom = [geom]
    end;

    # columns coordinates
    coord = self.col.coord

    if isa(color,String)
        if color === "cos2" 
            color_vec = vec(sum(Matrix{Float64}(self.col.cos2[!,2:end])[:,axis],dims=2)) 
            if legend_title === nothing legend_title = "Cos2" end;
        end;

        if color === "contrib" 
            color_vec = vec(sum(Matrix{Float64}(self.col.contrib[!,2:end])[:,axis],dims=2))
            if legend_title === nothing legend_title = "Contrib" end;
        end;

        if color in names(coord)[2:end]
            color_vec = vec(coord[!,color])
            if !isa(color_vec,Array{<:Number,1}) throw(ArgumentError("'color' must be a numeric variable.")) end;
            if legend_title === nothing legend_title = color end;
        end;
    elseif isa(color,Array{<:Number,1})
        color_vec = vec(color)
        if length(color_vec) !== nrow(coord) throw(ArgumentError("'color' must be a numeric vector of length $(nrow(coord)).")) end;
        if legend_title === nothing legend_title = "Cont_var" end;
    end;

    # Set x and y
    x_dim, y_dim = axis[1]+1, axis[2]+1

    # Set x labels
    if x_label === nothing x_label = "Dim." * string(axis[1]) *" ("*string(round.(self.eig[axis[1],4],digits=2))*"%)" end;
    # Set y labels
    if y_label === nothing y_label = "Dim." * string(axis[2]) *" ("*string(round.(self.eig[axis[2],4],digits=2))*"%)" end;
    # Set title
    if title === nothing title = "Columns Factor Map - CA" end;

    fig = Figure(size=figsize)
    ax = Axis(fig[1,1],xlabel=x_label,ylabel=y_label,title=title)

    if (isa(color,String) && (color in ["cos2","contrib"])) || isa(color,Array{<:Number,1})
        color_range = [minimum(color_vec),maximum(color_vec)]
        if :point in geom scatter!(ax,coord[!,x_dim],coord[!,y_dim],marker=marker,color=color_vec,colormap=palette,colorrange = color_range,markersize=point_size) end;
        if :text in geom text!(ax,coord[!,x_dim],coord[!,y_dim],text=coord[!,1],color=color_vec,colormap=palette,colorrange = color_range,fontsize=text_size,position=(ha,va)) end;
        if add_color_bar Colorbar(fig[1, 2],colorrange = color_range,colormap = palette, label = legend_title) end;
    else
        if :point in geom scatter!(ax,coord[!,x_dim],coord[!,y_dim],marker=marker,color=color,markersize=point_size) end;
        if :text in geom text!(ax,coord[!,x_dim],coord[!,y_dim],text=coord[!,1],color=color,fontsize=text_size,position=(ha,va)) end;
    end;
    
    # Add supplementary columns coordinates
    if col_sup
        if :col_sup in keys(self) 
            col_sup_coord = self.col_sup.coord
            if :point in geom scatter!(ax,col_sup_coord[!,x_dim],col_sup_coord[!,y_dim],marker=col_sup_marker,color=col_sup_color,markersize=col_sup_point_size)end;
            if :text in geom text!(ax,col_sup_coord[!,x_dim],col_sup_coord[!,y_dim],text=col_sup_coord[!,1],color=col_sup_color,fontsize=col_sup_text_size,position=(ha,va)) end;
        end;
    end;

    # Set x limits
    if x_lim !== nothing xlims!(ax,low = x_lim[1],high=x_lim[2]) end;
    # Set y limits
    if y_lim !== nothing ylims!(ax,low = y_lim[1],high=y_lim[2]) end;

    # Add horizontal line (hline)
    if add_hline hlines!(ax,0,color = hline_color,linestyle = hline_style,linewidth=hline_width) end;
    # Add vertical line (vline)
    if add_vline vlines!(ax,0,color = vline_color,linestyle = vline_style,linewidth=vline_width) end;
    # Hide spine
    if hide_spine hidespines!(ax) end;

    return fig
end;
"""
    fviz_ca_biplot(self;...)

Visualize Correspondence Analysis (CA) - Biplot of rows and columns

...
# Description
Correspondence analysis (CA) is an extension of Principal Component Analysis (PCA) suited to analyze frequencies formed by two categorical variables. `fviz_ca_biplot` provides CairoMakie based elegant visualization of CA outputs for rows and columns.

# Required arguments
- `self`: a CA NamedTuple

# Optional arguments
see `fviz_ca_row`, `fviz_ca_col`

# Return
a CairoMakie plot

# Examples
```julia-repl
julia> using Scientisttools
julia> children = get_dataset("children");
julia> res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
julia> # Biplot of rows and columns
julia> fig = fviz_ca_biplot(res_ca);
julia> fig
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fviz_ca_biplot(self::NamedTuple;
                         axis::AbstractVector=[1,2],
                         x_label=nothing,
                         y_label=nothing,
                         title="CA - Biplot",
                         x_lim=nothing,
                         y_lim=nothing,
                         row_marker::Symbol=:circle,
                         col_marker::Symbol=:diamond,
                         row_text_size=15,
                         col_text_size=15,
                         row_point_size=10,
                         col_point_size=10,
                         row_geom=[:point,:text],
                         col_geom=[:point,:text],
                         row_color::Symbol=:black,
                         col_color::Symbol=:steelblue,
                         row_sup::Bool=true,
                         row_sup_color::Symbol=:blue,
                         row_sup_marker::Symbol=:rect,
                         row_sup_text_size=15,
                         row_sup_point_size=10,
                         quali_sup::Bool=true,
                         quali_sup_color::Symbol=:red,
                         quali_sup_marker::Symbol=:diamond,
                         quali_sup_text_size=15,
                         quali_sup_point_size=10,
                         col_sup::Bool=true,
                         col_sup_color::Symbol=:green,
                         col_sup_marker::Symbol=:hexagon,
                         col_sup_text_size=15,
                         col_sup_point_size=10,
                         add_hline::Bool=true,
                         hline_color::Symbol=:black,
                         hline_style::Symbol=:dash,
                         hline_width=2,
                         add_vline::Bool=true,
                         vline_color::Symbol=:black,
                         vline_style::Symbol=:dash,
                         vline_width=2,
                         ha::Symbol=:center,
                         va::Symbol=:top,
                         hide_spine::Bool=false,
                         figsize=(800,500))

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if CA model
    if self.model !== "ca" throw(ArgumentError("self must be a CA NamedTuple.")) end;
    # Check axis
    if (length(axis) !== 2) || (axis[1] < 1) || (axis[2] > self.call.n_components) || (axis[1] > axis[2]) throw(ArgumentError("You must pass a valid 'axis'.")) end;

    # Extract rows and columnns coordinates
    row_coord, col_coord = self.row.coord, self.col.coord

    # Set x and y
    x_dim, y_dim = axis[1]+1, axis[2]+1

     # Set x labels
     if x_label === nothing x_label = "Dim." * string(axis[1]) *" ("*string(round.(self.eig[axis[1],4],digits=2))*"%)" end;
    # Set y labels
    if y_label === nothing y_label = "Dim." * string(axis[2]) *" ("*string(round.(self.eig[axis[2],4],digits=2))*"%)" end;
    # Set title
    if title === nothing title = "CA - Biplot" end;

    fig = Figure(size=figsize)
    ax = Axis(fig[1,1],xlabel=x_label,ylabel=y_label,title=title)

    #################################################################################################################################
    ## Rows informations
    #################################################################################################################################

    if :point in row_geom scatter!(ax,row_coord[!,x_dim],row_coord[!,y_dim],marker=row_marker,color=row_color,markersize=row_point_size) end;
    if :text in row_geom text!(ax,row_coord[!,x_dim],row_coord[!,y_dim],text=row_coord[!,1],color=row_color,fontsize=row_text_size,position=(ha,va)) end;
    
    # Add supplementary rows coordinates
    if row_sup
        if :row_sup in keys(self) 
            row_sup_coord = self.row_sup.coord
            if :point in row_geom scatter!(ax,row_sup_coord[!,x_dim],row_sup_coord[!,y_dim],marker=row_sup_marker,color=row_sup_color,markersize=row_sup_point_size)end;
            if :text in row_geom text!(ax,row_sup_coord[!,x_dim],row_sup_coord[!,y_dim],text=row_sup_coord[!,1],color=row_sup_color,fontsize=row_sup_text_size,position=(ha,va)) end;
        end;
    end;

    # Add supplementary qualitatives variables
    if quali_sup
        if :quali_sup in keys(self)
            quali_sup_coord = self.quali_sup.coord
            if :point in row_geom scatter!(ax,quali_sup_coord[!,x_dim],quali_sup_coord[!,y_dim],marker=quali_sup_marker,color=quali_sup_color,markersize=quali_sup_point_size)end;
            if :text in row_geom text!(ax,quali_sup_coord[!,x_dim],quali_sup_coord[!,y_dim],text=quali_sup_coord[!,1],color=quali_sup_color,fontsize=quali_sup_text_size,position=(ha,va)) end;
        end;
    end;

    #################################################################################################################################
    ## Columns informations
    #################################################################################################################################

    if :point in col_geom scatter!(ax,col_coord[!,x_dim],col_coord[!,y_dim],marker=col_marker,color=col_color,markersize=col_point_size) end;
    if :text in col_geom text!(ax,col_coord[!,x_dim],col_coord[!,y_dim],text=col_coord[!,1],color=col_color,fontsize=col_text_size,position=(ha,va)) end;
    
    # Add supplementary columns coordinates
    if col_sup
        if :col_sup in keys(self) 
            col_sup_coord = self.col_sup.coord
            if :point in col_geom scatter!(ax,col_sup_coord[!,x_dim],col_sup_coord[!,y_dim],marker=col_sup_marker,color=col_sup_color,markersize=col_sup_point_size)end;
            if :text in col_geom text!(ax,col_sup_coord[!,x_dim],col_sup_coord[!,y_dim],text=col_sup_coord[!,1],color=col_sup_color,fontsize=col_sup_text_size,position=(ha,va)) end;
        end;
    end;

    #################################################################################################################################
    ## Others informations
    #################################################################################################################################

    # Set x limits
    if x_lim !== nothing xlims!(ax,low = x_lim[1],high=x_lim[2]) end;
    # Set y limits
    if y_lim !== nothing ylims!(ax,low = y_lim[1],high=y_lim[2]) end;
    # Add horizontal line (hline)
    if add_hline hlines!(ax,0,color = hline_color,linestyle = hline_style,linewidth=hline_width) end;
    # Add vertical line (vline)
    if add_vline vlines!(ax,0,color = vline_color,linestyle = vline_style,linewidth=vline_width) end;
    # Hide spine
    if hide_spine hidespines!(ax) end;

    return fig
end;
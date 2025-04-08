"""
    fviz_pca_ind(self;...)

Visualize Principal Component Analysis (PCA) - Graph of individuals

...
# Description
Principal components analysis (PCA) reduces the dimensionality of multivariate data, to two or three that can be visualized graphically with minimal loss of information. `fviz_pca_ind` provides CairoMakie based elegant visualization of PCA outputs for individuals.

# Required arguments
- `self`: a PCA NamedTuple

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
- `ìnd_sup`: a boolean to either add or not supplementary individuals (by default = true).
- `ind_sup_color`: a symbol specifying a color for the supplementary individuals points (by default = :blue).
- `ind_sup_marker`: a marker style for the supplementary individuals points (by default = :rect).
- `ind_sup_point_size`: a numeric value specifying supplementary individuals marker size (by default = 10).
- `ind_sup_text_size`: a numeric value specifying supplementary individuals label size (by default = 15).
- `quali_sup`: a boolean to either add or not supplementary categories (by default = true).
- `quali_sup_color`: a symbol specifying a color for the supplementary categories points (by default = :red).
- `quali_sup_marker`: a marker style for the supplementary categories points (by default = :hexagon).
- `quali_sup_point_size`: a numeric value specifying supplementary categories marker size (by default = 10).
- `quali_sup_text_size`: a numeric value specifying supplementary categories label size (by default = 15).
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
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Graph of individuals
julia> fig = fviz_pca_ind(res_pca);
julia> fig
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fviz_pca_ind(self::NamedTuple;
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
                      ind_sup::Bool=true,
                      ind_sup_color::Symbol=:blue,
                      ind_sup_marker::Symbol=:rect,
                      ind_sup_point_size=10,
                      ind_sup_text_size=15,
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
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;
    # Check axis
    if (length(axis) !== 2) || (axis[1] < 1) || (axis[2] > self.call.n_components) || (axis[1] > axis[2])  throw(ArgumentError("You must pass a valid 'axis'.")) end;

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
    coord = hcat(self.ind.coord,X)

    # Add supplementary quantitative variables
    if :quanti_sup in keys(self)
        X_quanti_sup = self.call.Xtot[!,self.call.quanti_sup]
        if isa(X_quanti_sup,AbstractVector) X_quanti_sup = DataFrame(self.call.quanti_sup => X_quanti_sup) end;
        if self.call.ind_sup !== nothing X_quanti_sup = X_quanti_sup[Not(self.call.ind_sup),:] end;
        X_quanti_sup = convert_to_float(X_quanti_sup)
        coord = hcat(coord,X_quanti_sup)
    end;

    # Add supplementary quantitative variables
    if :quali_sup in keys(self)
        X_quali_sup = self.call.Xtot[!,self.call.quali_sup]
        if isa(X_quali_sup,AbstractVector) X_quali_sup = DataFrame(self.call.quali_sup => X_quali_sup) end;
        if self.call.ind_sup !== nothing X_quali_sup = X_quali_sup[Not(self.call.ind_sup),:] end;
        coord = hcat(coord,X_quali_sup)
    end;

    # For
    if isa(color,String)
        if color === "cos2" 
            color_vec = vec(sum(Matrix{Float64}(self.ind.cos2[!,2:end])[:,axis],dims=2)) 
            if legend_title === nothing legend_title = "Cos2" end;
        end;

        if color === "contrib" 
            color_vec = vec(sum(Matrix{Float64}(self.ind.contrib[!,2:end])[:,axis],dims=2))
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
    if title === nothing title = "Individuals Factor Map - PCA" end;

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
        axislegend(ax,legend_title)
    end;

    # Add supplementary individuals coordinates
    if ind_sup
        if :ind_sup in keys(self) 
            ind_sup_coord = self.ind_sup.coord
            if :point in geom scatter!(ax,ind_sup_coord[!,x_dim],ind_sup_coord[!,y_dim],marker=ind_sup_marker,color=ind_sup_color,markersize=ind_sup_point_size)end;
            if :text in geom text!(ax,ind_sup_coord[!,x_dim],ind_sup_coord[!,y_dim],text=ind_sup_coord[!,1],color=ind_sup_color,fontsize=ind_sup_text_size,position=(ha,va)) end;
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
    fviz_pca_var(self;...)

Visualize Principal Component Analysis (PCA) - Graph of variables

...
# Description
Principal components analysis (PCA) reduces the dimensionality of multivariate data, to two or three that can be visualized graphically with minimal loss of information. `fviz_pca_var` provides CairoMakie based elegant visualization of PCA outputs for variables.

# Required arguments
- `self`: a PCA NamedTuple

# Optional arguments
- `axis`: a numeric vector/tuple of length 2 specifying the dimension to be plotted (by default =[0,1]).
- `x_lim`: a numeric vector/tuple of length 2 specifying the range of plotted 'x' values (by default = (-1.2,1.2)).
- `y_lim`: a numeric vector/tuple of length 2 specifying the range of plotted 'y' values (by default = (-1.2,1.2)).
- `x_label`: a string specifying the label text of x (by default = nothing).
- `y_label`: a string specifying the label text of y (by default = nothing).
- `title`: a string corresponding to the title of the graph you draw (by default = nothing).
- `color`: a symbol specifying the color for active variables (by default = :black)
- `geom`: a symbol specifying the geometry to be used for the graph. Allowed values are the combinaison of [:arrow,:text]. Use ':arrow'  (to show only arrows); ':text' to show only labels; [:arrow,:text] to show both types.
- `arrow_width`: a numeric specifying the width of the arrow (by default = 1)
- `arrow_size`: a numeric specifying the size of the arrow (by default = 10)
- `line_style`: a symbol specifying the linestyle of the arrow of active variables (by default = :solid)
- `palette`: a symbol specifying the palette color of the colorbar (by default = :Hiroshige).
- `text_size`: a numeric value specifying the label size (by default = 15).
- `add_color_bar`: a boolean, default = true. If `true`, the colorbar is plotted, else the colorbar isn't plotted.
- `legend_title`: a string corresponding to the title of the legend (by default = nothing).
- `habillage`: a string corresponding the name of group variable (by default = nothing).
- `quanti_sup`: a boolean to either add or not supplementary quantitative variables (by default = true).
- `quanti_sup_color`: a symbol specifying a color for the supplementary quuantitative variables (by default = :blue).
- `quanti_sup_arrow_width`: a numeric specifying the width of the supplementary arrow (by default = 1)
- `quanti_sup_arrow_size`: a numeric specifying the size of the supplementary arrow (by default = 10)
- `quanti_sup_line_style`: a symbol specifying the linestyle of the arrow of supplementary quantitative variables (by default = :dash)
- `quanti_sup_text_size`: a numeric value specifying the supplementary label size (by default = 15).
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
- `add_circle`: a boolean, whether to add or not a circle to plot (by default = true).
- `color_circle` : a symbol specifying the color for the correlation circle (by default = :gray)
- `figsize`: the aspect of figure (by default = (500,500)).

# Return
a CairoMakie plot

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Graph of variables
julia> fig = fviz_pca_var(res_pca);
julia> fig
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fviz_pca_var(self::NamedTuple;
                      axis::AbstractVector=[1,2],
                      x_lim=(-1.2,1.2),
                      y_lim=(-1.2,1.2),
                      x_label=nothing,
                      y_label=nothing,
                      title = nothing,
                      color=:black,
                      geom=[:arrow,:text],
                      arrow_width=1,
                      arrow_size=10,
                      line_style::Symbol=:solid,
                      palette::Symbol=:Hiroshige,
                      text_size=15,
                      add_color_bar::Bool=true,
                      legend_title=nothing,
                      habillage=nothing,
                      quanti_sup::Bool=true,
                      quanti_sup_arrow_width=1,
                      quanti_sup_arrow_size=10,
                      quanti_sup_color::Symbol=:blue,
                      quanti_sup_line_style::Symbol=:dash,
                      quanti_sup_text_size=15,
                      add_hline::Bool=true,
                      hline_color::Symbol=:black,
                      hline_style::Symbol=:dash,
                      hline_width=2,
                      add_vline::Bool=true,
                      vline_color::Symbol=:black,
                      vline_style::Symbol=:dash,
                      vline_width=2,
                      ha::Symbol=:left,
                      va::Symbol=:top,
                      hide_spine::Bool=false, 
                      add_circle::Bool=true,
                      color_circle::Symbol=:black,
                      figsize=(500,500))

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;
    # Check axis
    if (length(axis) !== 2) || (axis[1] < 1) || (axis[2] > self.call.n_components) || (axis[1] > axis[2]) throw(ArgumentError("You must pass a valid 'axis'.")) end;

    #
    if isa(geom,AbstractVector)
        intersect = [x for x in geom if x in [:arrow,:text]]
        if length(intersect) == 0 throw(ArgumentError("Specified value(s) for the argument geom aren't allowed.")) end;
    elseif Base.isidentifier(geom)
        if !(geom in [:point,:text]) throw(ArgumentError("Invalid geom_type name $(geom)")) end;
        geom = [geom]
    end;

    # Add active data
    coord = self.var.coord

    # Set clor
    if isa(color,String)
        if color === "cos2" 
            color_vec = vec(sum(Matrix{Float64}(self.var.cos2[!,2:end])[:,axis],dims=2)) 
            if legend_title === nothing legend_title = "Cos2" end;
        end;

        if color === "contrib" 
            color_vec = vec(sum(Matrix{Float64}(self.var.contrib[!,2:end])[:,axis],dims=2))
            if legend_title === nothing legend_title = "Contrib" end;
        end;
    elseif isa(color,Array{<:Number,1})
        color_vec = vec(color)
        if length(color_vec) !== nrow(coord) throw(ArgumentError("'color' must be a numeric vector of length $(nrow(coord)).")) end;
        if legend_title === nothing legend_title = "Cont_var" end;
    end;

    # For habillage
    if (habillage !== nothing) && isa(habillage,Array{String,1})
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
    if title === nothing title = "Variables Factor Map - PCA" end;

    fig = Figure(size=figsize)
    ax = Axis(fig[1,1],xlabel=x_label,ylabel=y_label,title=title)

    if habillage === nothing
        if (isa(color,String) && (color in ["cos2","contrib"])) || isa(color,Array{<:Number,1})
            color_range = [minimum(color_vec),maximum(color_vec)]
            if :arrow in geom arrows!(ax,zeros(nrow(coord)),zeros(nrow(coord)),coord[!,x_dim], coord[!,y_dim],
                                    color=color_vec,arrowcolor=color_vec,colormap=palette,colorrange = color_range,
                                    linewidth=arrow_width,arrowsize=arrow_size,linestyle=line_style) end;
            if :text in geom text!(ax,coord[!,x_dim],coord[!,y_dim],text=coord[!,1],color=color_vec,colormap=palette,
                                colorrange = color_range,fontsize=text_size,position=(ha,va)) end;
            if add_color_bar Colorbar(fig[1, 2],colorrange = color_range,colormap = palette, label = legend_title) end;
        else
            if :arrow in geom arrows!(ax,zeros(nrow(coord)),zeros(nrow(coord)), coord[!,x_dim],coord[!,y_dim],color=color,
                                    arrowcolor=color,linewidth=arrow_width,arrowsize=arrow_size,linestyle=line_style) end;
            if :text in geom text!(ax,coord[!,x_dim],coord[!,y_dim],text=coord[!,1],color=color,fontsize=text_size,position=(ha,va)) end;
        end;
    else
        uq_group,colors = sort(unique(vsqual)),[:blue,:red,:green,:steelblue]
        Random.seed!(1234);
        color_group = Dict(zip(uq_group,sample(colors, length(uq_group); replace=false)))
        for g in uq_group
            idx = [i for i in 1:length(vsqual) if vsqual[i]==g]
            if :arrow in geom arrows!(ax,zeros(length(idx)),zeros(length(idx)),coord[idx,x_dim],coord[idx,y_dim],color=color_group[g],label=g,
                                      arrowcolor=color_group[g],linewidth=arrow_width,arrowsize=arrow_size,linestyle=line_style) end;
            if :text in geom text!(ax,coord[idx,x_dim],coord[idx,y_dim],text=coord[idx,1],color=color_group[g],fontsize=text_size,position=(ha,va)) end;
        end;
        axislegend(ax,legend_title)
    end;

    # Add supplementary quantitative variable coordinates
    if quanti_sup
        if :quanti_sup in keys(self)
            quanti_sup_coord = self.quanti_sup.coord
            if :arrow in geom arrows!(ax,zeros(nrow(quanti_sup_coord)),zeros(nrow(quanti_sup_coord)),quanti_sup_coord[!,x_dim],quanti_sup_coord[!,y_dim],
                                      color=quanti_sup_color,arrowcolor=quanti_sup_color,linewidth=quanti_sup_arrow_width,arrowsize=quanti_sup_arrow_size,linestyle=quanti_sup_line_style) end;
            if :text in geom text!(ax,quanti_sup_coord[!,x_dim],quanti_sup_coord[!,y_dim],text=quanti_sup_coord[!,1],color=quanti_sup_color,fontsize=quanti_sup_text_size,position=(ha,va)) end;
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
    # Add circle
    if add_circle arc!(ax,Point2f(0), 1, -π, π,color=color_circle) end;
    # Hide spine
    if hide_spine hidespines!(ax) end;

    return fig
end;
"""
    fviz_pca_biplot(self;...)

Visualize Principal Component Analysis (PCA) - Biplot of individuals and variables

...
# Description
Principal components analysis (PCA) reduces the dimensionality of multivariate data, to two or three that can be visualized graphically with minimal loss of information. `fviz_pca_biplot` provides CairoMakie based elegant visualization of PCA outputs for individuals and variables.

# Required arguments
- `self`: a PCA NamedTuple

# Optional arguments
see `fviz_pca_ind`, `fviz_pca_var`

# Return
a CairoMakie plot

# Examples
```julia-repl
julia> using Scientisttools
julia> decathlon = get_dataset("decathlon");
julia> res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
julia> # Biplot of individuals and variables
julia> fig = fviz_pca_biplot(res_pca);
julia> fig
```

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function fviz_pca_biplot(self::NamedTuple;
                         axis::AbstractVector=[1,2],
                         x_label=nothing,
                         y_label=nothing,
                         title="PCA - Biplot",
                         x_lim=nothing,
                         y_lim=nothing,
                         marker::Symbol=:circle,
                         ind_text_size=15,
                         var_text_size=15,
                         ind_point_size=10,
                         ind_geom=[:point,:text],
                         var_geom=[:arrow,:text],
                         ind_color::Symbol=:black,
                         var_color::Symbol=:steelblue,
                         arrow_width=1,
                         arrow_size=10,
                         line_style::Symbol=:solid,
                         add_circle::Bool=false,
                         color_circle::Symbol=:gray,
                         ind_sup::Bool=true,
                         ind_sup_color::Symbol=:blue,
                         ind_sup_marker::Symbol=:rect,
                         ind_sup_text_size=15,
                         ind_sup_point_size=10,
                         quanti_sup::Bool=true,
                         quanti_sup_color::Symbol=:green,
                         quanti_sup_text_size=15,
                         quanti_sup_arrow_width=1,
                         quanti_sup_arrow_size=10,
                         quanti_sup_line_style::Symbol=:dash,
                         quali_sup::Bool=true,
                         quali_sup_color::Symbol=:red,
                         quali_sup_marker::Symbol=:diamond,
                         quali_sup_text_size=15,
                         quali_sup_point_size=10,
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
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;
    # Check axis
    if (length(axis) !== 2) || (axis[1] < 1) || (axis[2] > self.call.n_components) || (axis[1] > axis[2])  throw(ArgumentError("You must pass a valid 'axis'.")) end;

    # Individuals coordinates
    ind_names, ind_coord = self.ind.coord[!,[1]], Matrix(self.ind.coord[!,2:end])[:,axis]
    var_names, var_coord = self.var.coord[!,[1]], Matrix(self.var.coord[!,2:end])[:,axis]

    # Rescale variables coordinates
    xscale = (maximum(ind_coord[:,1]) - minimum(ind_coord[:,1]))/(maximum(var_coord[:,1]) - minimum(var_coord[:,1]))
    yscale = (maximum(ind_coord[:,2]) - minimum(ind_coord[:,2]))/(maximum(var_coord[:,2]) - minimum(var_coord[:,2]))
    rscale = min(xscale,yscale)

    #
    ind_coord = hcat(ind_names,DataFrame(ind_coord,self.call.dim_index[axis]))
    var_coord = hcat(var_names,DataFrame(var_coord * rscale, self.call.dim_index[axis]))

    # Set x and y
    x_dim, y_dim = axis[1]+1, axis[2]+1

     # Set x labels
     if x_label === nothing x_label = "Dim." * string(axis[1]) *" ("*string(round.(self.eig[axis[1],4],digits=2))*"%)" end;
    # Set y labels
    if y_label === nothing y_label = "Dim." * string(axis[2]) *" ("*string(round.(self.eig[axis[2],4],digits=2))*"%)" end;
    # Set title
    if title === nothing title = "PCA - Biplot" end;

    fig = Figure(size=figsize)
    ax = Axis(fig[1,1],xlabel=x_label,ylabel=y_label,title=title)

    #################################################################################################################################
    ## Individuals informations
    #################################################################################################################################

    if :point in ind_geom scatter!(ax,ind_coord[!,x_dim],ind_coord[!,y_dim],marker=marker,color=ind_color,markersize=ind_point_size) end;
    if :text in ind_geom text!(ax,ind_coord[!,x_dim],ind_coord[!,y_dim],text=ind_coord[!,1],color=ind_color,fontsize=ind_text_size,position=(ha,va)) end;
    
    # Add supplementary individuals coordinates
    if ind_sup
        if :ind_sup in keys(self) 
            ind_sup_coord = self.ind_sup.coord
            if :point in ind_geom scatter!(ax,ind_sup_coord[!,x_dim],ind_sup_coord[!,y_dim],marker=ind_sup_marker,color=ind_sup_color,markersize=ind_sup_point_size)end;
            if :text in ind_geom text!(ax,ind_sup_coord[!,x_dim],ind_sup_coord[!,y_dim],text=ind_sup_coord[!,1],color=ind_sup_color,fontsize=ind_sup_text_size,position=(ha,va)) end;
        end;
    end;

    # Add supplementary qualitatives variables
    if quali_sup
        if :quali_sup in keys(self)
            quali_sup_coord = self.quali_sup.coord
            if :point in ind_geom scatter!(ax,quali_sup_coord[!,x_dim],quali_sup_coord[!,y_dim],marker=quali_sup_marker,color=quali_sup_color,markersize=quali_sup_point_size)end;
            if :text in ind_geom text!(ax,quali_sup_coord[!,x_dim],quali_sup_coord[!,y_dim],text=quali_sup_coord[!,1],color=quali_sup_color,fontsize=quali_sup_text_size,position=(ha,va)) end;
        end;
    end;

    #################################################################################################################################
    ## Variables informations
    #################################################################################################################################

    if :arrow in var_geom arrows!(ax,zeros(nrow(var_coord)),zeros(nrow(var_coord)), var_coord[!,x_dim],var_coord[!,y_dim],color=var_color,
                                  arrowcolor=var_color,linewidth=arrow_width,arrowsize=arrow_size,linestyle=line_style) end;
    if :text in var_geom text!(ax,var_coord[!,x_dim],var_coord[!,y_dim],text=var_coord[!,1],color=var_color,fontsize=var_text_size,position=(ha,va)) end;
    
    # Add supplementary quantitative variable coordinates
    if quanti_sup
        if :quanti_sup in keys(self)
            quanti_sup_names = self.quanti_sup.coord[!,[1]]
            quanti_sup_coord = hcat(quanti_sup_names,DataFrame(Matrix(self.quanti_sup.coord[!,2:end])[:,axis] * rscale,self.call.dim_index[axis]))
            if :arrow in var_geom arrows!(ax,zeros(nrow(quanti_sup_coord)),zeros(nrow(quanti_sup_coord)),quanti_sup_coord[!,x_dim],quanti_sup_coord[!,y_dim],
                                          color=quanti_sup_color,arrowcolor=quanti_sup_color,linewidth=quanti_sup_arrow_width,arrowsize=quanti_sup_arrow_size,linestyle=quanti_sup_line_style) end;
            if :text in var_geom text!(ax,quanti_sup_coord[!,x_dim],quanti_sup_coord[!,y_dim],text=quanti_sup_coord[!,1],color=quanti_sup_color,fontsize=quanti_sup_text_size,position=(ha,va)) end;
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
    # Add circle
    if add_circle arc!(ax,Point2f(0), 1, -π, π,color=color_circle) end;
    # Hide spine
    if hide_spine hidespines!(ax) end;

    return fig
end;
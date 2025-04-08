"""
    fviz_pca_ind(self;...)

Visualize Principal Component Analysis (PCA) - Graph of individuals

...
# Description
Principal components analysis (PCA) reduces the dimensionality of multivariate data, to two or three that can be visualized graphically with minimal loss of information. `fviz_pca_ind` provides TidierPlots based elegant visualization of PCA outputs for individuals.

# Required arguments
- `self`: a PCA NamedTuple

# Optional arguments
- `axis`: a vector of 2 
...
"""
function fviz_pca_ind(self::NamedTuple;
                      axis::AbstractVector=[1,2],
                      x_lim=nothing,
                      y_lim=nothing,
                      x_label=nothing,
                      y_label=nothing,
                      title = nothing,
                      color = "black",
                      palette = "Hiroshige",
                      habillage = nothing,
                      geom=[:point,:text],
                      ind_sup::Bool=true,
                      color_sup::String="blue",
                      quali_sup::Bool=true,
                      color_quali_sup::String="red",
                      add_hline::Bool=true,
                      add_vline::Bool=true,
                      ggtheme=nothing)

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;

    if length(axis) !== 2 || axis[1] < 1 || axis[2] > self.call.n_components || axis[1] > axis[2] 
        throw(ArgumentError("You must pass a valid 'axis'.")) 
    end;

    #
    if isa(geom,AbstractVector)
        intersect = [x for x in geom if x in [:point,:text]]
        if length(intersect) == 0 throw(ArgumentError("Specified value(s) for the argument geom aren't allowed.")) end;
    elseif Base.isidentifier(geom)
        if !(geom in [:point,:text]) throw(ArgumentError("Invalid geom_type name $(geom)")) end;
        geom = [geom]
    end;

    # Add active data
    coord = leftjoin(self.ind.coord,self.call.X,on=:Individuals)

    # Add supplementary quantitative variables
    if :quanti_sup in keys(self)
        X_quanti_sup = self.call.Xtot[!,self.call.quanti_sup]
        if isa(X_quanti_sup,AbstractVector) X_quanti_sup = DataFrame(self.call.quanti_sup => X_quanti_sup) end;
        if self.call.ind_sup !== nothing X_quanti_sup = X_quanti_sup[Not(self.call.ind_sup),:] end;
        coord = hcat(coord,X_quanti_sup)
    end;

    # Add supplementary quantitative variables
    if :quali_sup in keys(self)
        X_quali_sup = self.call.Xtot[!,self.call.quali_sup]
        if isa(X_quali_sup,AbstractVector) X_quali_sup = DataFrame(self.call.quali_sup => X_quali_sup) end;
        if self.call.ind_sup !== nothing X_quali_sup = X_quali_sup[Not(self.call.ind_sup),:] end;
        coord = hcat(coord,X_quali_sup)
    end;

    #
    if isa(color,String) && (color in ["cos2","contrib"])
        if color === "cos2" 
            color_data = DataFrame("cos2" => vec(sum(Matrix{Float64}(self.ind.cos2[!,2:end])[:,axis],dims=2))) 
            coord = hcat(coord,color_data)
            color_name = "Cos2"
        end;

        if color === "contrib" 
            color_data = DataFrame("contrib" => vec(sum(Matrix{Float64}(self.ind.contrib[!,2:end])[:,axis],dims=2))) 
            coord = hcat(coord,color_data)
            color_name = "Contrib"
        end;
    end;

    # Set x and y
    x_dim, y_dim = "Dim." * string(axis[1]), "Dim." * string(axis[2])

    # Initialisation
    fig = gg.ggplot(coord)
    if habillage === nothing
        if (isa(color,String) && (color in reduce(vcat,[["cos2","contrib"],names(coord)[2:end]])))
            if :point in geom fig = fig + gg.geom_point(gg.aes(x = x_dim, y = y_dim, color=color)) end;
            if :text in geom fig = fig + gg.geom_text(gg.aes(x = x_dim,y = y_dim,label="Individuals",color=color)) end;
            fig = fig + gg.scale_color_continuous(palette = palette)
        else
            if :point in geom fig = fig + gg.geom_point(gg.aes(x = x_dim, y = y_dim), color = color) end;
            if :text in geom fig = fig + gg.geom_text(gg.aes(x = x_dim, y = y_dim, label="Individuals"), color=color) end;
        end;
    else
        if !(habillage in names(coord)[2:end]) throw(ArgumentError("$habillage not in DataFrame.")) end;
        if :point in geom fig = fig + gg.geom_point(gg.aes(x = x_dim, y = y_dim,color=habillage)) end;
        if :text in geom fig = fig + gg.geom_text(gg.aes(x = x_dim, y = y_dim, label="Individuals",color=habillage)) end;
    end;

    # Add supplementary individuals coordinates
    if ind_sup
        if :ind_sup in keys(self) 
            ind_sup_coord = self.ind_sup.coord
            if :point in geom fig = fig + gg.geom_point(data=ind_sup_coord,gg.aes(x = x_dim, y = y_dim), color = color_sup) end;
            if :text in geom fig = fig + gg.geom_text(data=ind_sup_coord,gg.aes(x = x_dim, y = y_dim, label="Individuals"), color=color_sup) end;
        end;
    end;

    # Add supplementary qualitatives variables
    if quali_sup
        if :quali_sup in keys(self)
            if habillage === nothing
                quali_sup_coord = self.quali_sup.coord
                if :point in geom fig = fig + gg.geom_point(data=quali_sup_coord,gg.aes(x = x_dim, y = y_dim), color = color_quali_sup) end;
                if :text in geom fig = fig + gg.geom_text(data=quali_sup_coord,gg.aes(x = x_dim, y = y_dim, label="Modalities"), color=color_quali_sup) end;
            end;
        end;
    end;

    # Set x labels
    if x_label === nothing x_label = x_dim *" ("*string(round.(self.eig[axis[1],4],digits=2))*"%)" end;
    # Set y labels
    if y_label === nothing y_label = y_dim *" ("*string(round.(self.eig[axis[2],4],digits=2))*"%)" end;
    # Set title
    if title === nothing title = "Individuals Factor Map - PCA" end;
    fig = fig + gg.labs(x=x_label, y=y_label, title=title)

    # Set x limits
    if x_lim === nothing x_lim = (minimum(coord[!,axis[1]+1]),maximum(coord[!,axis[1]+1])) end;
    # Set y limits
    if y_lim === nothing y_lim = (minimum(coord[!,axis[2]+1]),maximum(coord[!,axis[2]+1])) end;
    fig = fig + gg.lims(x = x_lim, y = y_lim)

    # Add horizontal line (hline)
    if add_hline fig = fig + gg.geom_hline(yintercept=0) end;
    # Add vertical line (vline)
    if add_vline fig = fig + gg.geom_vline(xintercept=0) end;

    # Add theme
    if ggtheme === nothing ggtheme = gg.theme_ggplot2() end;
    return fig + ggtheme
end;

"""
    fviz_pca_var(self;...)

Visualize Principal Component Analysis (PCA) - Graph of variables

...
# Description
Principal components analysis (PCA) reduces the dimensionality of multivariate data, to two or three that can be visualized graphically with minimal loss of information. `fviz_pca_var` provides TidierPlots based elegant visualization of PCA outputs for variables.

...
"""
function fviz_pca_var(self::NamedTuple;
                      axis::AbstractVector=[1,2],
                      x_label=nothing,
                      y_label=nothing,
                      title = nothing,
                      color = "black",
                      palette = "Hiroshige",
                      geom=[:arrow,:text],
                      quanti_sup::Bool=true,
                      color_sup::String="blue",
                      add_hline::Bool=true,
                      add_vline::Bool=true,
                      add_circle::Bool=true,
                      color_circle="gray",
                      ggtheme=nothing)

    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;

    if length(axis) !== 2 || axis[1] < 1 || axis[2] > self.call.n_components || axis[1] > axis[2] 
        throw(ArgumentError("You must pass a valid 'axis'.")) 
    end;

    #
    if isa(geom,AbstractVector)
        intersect = [x for x in geom if x in [:arrow,:text]]
        if length(intersect) == 0 throw(ArgumentError("Specified value(s) for the argument geom aren't allowed.")) end;
    elseif Base.isidentifier(geom)
        if !(geom in [:arrow,:text]) throw(ArgumentError("Invalid geom_type name $(geom)")) end;
        geom = [geom]
    end;

    # Add active data
    coord = self.var.coord

    if (isa(color,String) && (color in ["cos2","contrib"]))
        if color === "cos2" 
            color_data = DataFrame("cos2" => vec(sum(Matrix{Float64}(self.var.cos2[!,2:end])[:,axis],dims=2))) 
            coord = hcat(coord,color_data)
            color_name = "Cos2"
        end;

        if color === "contrib" 
            color_data = DataFrame("contrib" => vec(sum(Matrix{Float64}(self.var.contrib[!,2:end])[:,axis],dims=2))) 
            coord = hcat(coord,color_data)
            color_name = "Contrib"
        end;
    end;

    # Set x and y
    x_dim, y_dim = "Dim." * string(axis[1]), "Dim." * string(axis[2])

    # Initialisation
    fig = gg.ggplot(coord)
    if (isa(color,String) && (color in ["cos2","contrib"]))
        #if :point in geom fig = fig + gg.geom_point(gg.aes(x = x_dim, y = y_dim, color=color)) end;
        if :text in geom fig = fig + gg.geom_text(gg.aes(x = x_dim,y = y_dim,label="Variables",color=color)) end;
        fig = fig + gg.scale_color_continuous(palette = palette)
    else
        if :arrow in geom fig = fig + geom_segment(gg,coord[!,[x_dim,y_dim]],color = color) end;
        if :text in geom fig = fig + gg.geom_text(gg.aes(x = x_dim, y = y_dim, label="Variables"), color=color) end;
    end;

    # Add supplementary quantitative variable coordinates
    if quanti_sup
        if :quanti_sup in keys(self) 
            quanti_sup_coord = self.quanti_sup.coord
            if :arrow in geom fig = fig + geom_segment(gg,quanti_sup_coord[!,[x_dim,y_dim]],color = color_sup) end;
            if :text in geom fig = fig + gg.geom_text(data=quanti_sup_coord,gg.aes(x = x_dim, y = y_dim, label="Variables"), color=color_sup) end;
        end;
    end;

    # Add circle
    if add_circle fig = fig + geom_circle(gg,x_dim,y_dim,color=color_circle) end;

    # Set x labels
    if x_label === nothing x_label = x_dim *" ("*string(round.(self.eig[axis[1],4],digits=2))*"%)" end;
    # Set y labels
    if y_label === nothing y_label = y_dim *" ("*string(round.(self.eig[axis[2],4],digits=2))*"%)" end;
    # Set title
    if title === nothing title = "Variables Factor Map - PCA" end;
    fig = fig + gg.labs(x=x_label, y=y_label, title=title)

    # Set x and y limits
    fig = fig + gg.lims(x = (-1.1,1.1), y = (-1.1,1.1))

    # Add horizontal line (hline)
    if add_hline fig = fig + gg.geom_hline(yintercept=0) end;
    # Add vertical line (vline)
    if add_vline fig = fig + gg.geom_vline(xintercept=0) end;

    # Add theme
    if ggtheme === nothing ggtheme = gg.theme_ggplot2() end;
    return fig + ggtheme
end;
"""
    fviz_pca_biplot(self;...)

Visualize Principal Component Analysis (PCA) - Biplot of individuals and variables

...
# Description
Principal components analysis (PCA) reduces the dimensionality of multivariate data, to two or three that can be visualized graphically with minimal loss of information. fviz_pca_biplot provides TidierPlots based elegant visualization of PCA outputs for individuals and variables.
...
"""
function fviz_pca_biplot(self::NamedTuple;
                         axis::AbstractVector=[1,2],
                         x_label=nothing,
                         y_label=nothing,
                         title="PCA - Biplot",
                         x_lim=nothing,
                         y_lim=nothing,
                         ind_geom=[:point,:text],
                         var_geom=[:arrow,:text],
                         ind_color="black",
                         var_color="steelblue",
                         add_circle::Bool=false,
                         color_circle="gray",
                         ind_sup::Bool=true,
                         ind_sup_color::String="blue",
                         quali_sup::Bool=true,
                         quali_sup_color::String="red",
                         quanti_sup::Bool=true,
                         quanti_sup_color::String="green",
                         add_hline::Bool=true,
                         add_vline::Bool=true,
                         ggtheme=nothing)
        
    # Check if self is an object of NamedTuple
    if !isa(self,NamedTuple) throw(ArgumentError("self must be a NamedTuple object.")) end;
    # Check if PCA model
    if self.model !== "pca" throw(ArgumentError("self must be a PCA NamedTuple.")) end;

    if length(axis) !== 2 || axis[1] < 1 || axis[2] > self.call.n_components || axis[1] > axis[2] 
        throw(ArgumentError("You must pass a valid 'axis'.")) 
    end;

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
    x_dim, y_dim = "Dim." * string(axis[1]), "Dim." * string(axis[2])

    # Initialize
    fig = gg.ggplot()

    #################################################################################################################################
    ## Individuals informations
    #################################################################################################################################

    if :point in ind_geom fig = fig + gg.geom_point(data=ind_coord,gg.aes(x = x_dim, y = y_dim), color = ind_color) end;
    if :text in ind_geom fig = fig + gg.geom_text(data=ind_coord,gg.aes(x = x_dim, y = y_dim, label="Individuals"), color=ind_color) end;

    # Add supplementary individuals coordinates
    if ind_sup
        if :ind_sup in keys(self) 
            ind_sup_coord = self.ind_sup.coord
            if :point in ind_geom fig = fig + gg.geom_point(data=ind_sup_coord,gg.aes(x = x_dim, y = y_dim), color = ind_sup_color) end;
            if :text in ind_geom fig = fig + gg.geom_text(data=ind_sup_coord,gg.aes(x = x_dim, y = y_dim, label="Individuals"), color=ind_sup_color) end;
        end;
    end;

    # Add supplementary qualitatives variables
    if quali_sup
        if :quali_sup in keys(self)
            quali_sup_coord = self.quali_sup.coord
            if :point in ind_geom fig = fig + gg.geom_point(data=quali_sup_coord,gg.aes(x = x_dim, y = y_dim), color = quali_sup_color) end;
            if :text in ind_geom fig = fig + gg.geom_text(data=quali_sup_coord,gg.aes(x = x_dim, y = y_dim, label="Modalities"), color=quali_sup_color) end;
        end;
    end;

    #################################################################################################################################
    ## Variables informations
    #################################################################################################################################

    if :arrow in var_geom fig = fig + geom_segment(gg,var_coord[!,[x_dim,y_dim]],color = var_color) end;
    if :text in var_geom fig = fig + gg.geom_text(data=var_coord,gg.aes(x = x_dim, y = y_dim, label="Variables"), color=var_color) end;

    # Add supplementary quantitative variable coordinates
    if quanti_sup
        if :quanti_sup in keys(self) 
            quanti_sup_names = self.quanti_sup.coord[!,[1]]
            quanti_sup_coord = hcat(quanti_sup_names,DataFrame(Matrix(self.quanti_sup.coord[!,2:end])[:,axis] * rscale,self.call.dim_index[axis]))
            if :arrow in var_geom fig = fig + geom_segment(gg,quanti_sup_coord[!,[x_dim,y_dim]],color = quanti_sup_color) end;
            if :text in var_geom fig = fig + gg.geom_text(data=quanti_sup_coord,gg.aes(x = x_dim, y = y_dim, label="Variables"), color=quanti_sup_color) end;
        end;
    end;

    # Add circle
    if add_circle fig = fig + geom_circle(gg,x_dim,y_dim,color=color_circle) end;

    # Set x labels
    if x_label === nothing x_label = x_dim *" ("*string(round.(self.eig[axis[1],4],digits=2))*"%)" end;
    # Set y labels
    if y_label === nothing y_label = y_dim *" ("*string(round.(self.eig[axis[2],4],digits=2))*"%)" end;
    # Set title
    if title === nothing title = "PCA - Biplot" end;
    fig = fig + gg.labs(x=x_label, y=y_label, title=title)

    # Set x limits
    if x_lim === nothing x_lim = (minimum(self.ind.coord[!,axis[1]+1]),maximum(self.ind.coord[!,axis[1]+1])) end;
    # Set y limits
    if y_lim === nothing y_lim = (minimum(self.ind.coord[!,axis[2]+1]),maximum(self.ind.coord[!,axis[2]+1])) end;
    fig = fig + gg.lims(x = x_lim, y = y_lim)

    # Add horizontal line (hline)
    if add_hline fig = fig + gg.geom_hline(yintercept=0) end;
    # Add vertical line (vline)
    if add_vline fig = fig + gg.geom_vline(xintercept=0) end;

    # Add theme
    if ggtheme === nothing ggtheme = gg.theme_ggplot2() end;
    return fig + ggtheme
end;
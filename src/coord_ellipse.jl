# Load packages
using DataFrames
using Distributions
using Statistics

function coord_ellipse(X::AbstractDataFrame;centre=nothing,axes=(1,2),level_conf::AbstractFloat=0.95,npoints::Int=100,bary::Bool=false)

    function ellipse(x;sigma=nothing,centre=(0,0),level::AbstractFloat=0.95,t = sqrt(quantile(Chisq(2),level)),which_index=(1,2),npoints::Int=100)
        # Exctract index
        if x isa AbstractMatrix
            x_idx, y_idx = which_index
            r = x[x_idx,y_idx]
            if sigma === nothing
                sigma = sqrt.((x[x_idx,x_idx],x[y_idx,y_idx]))
                if sigma[1] > 0 r = r/sigma[1] end;
                if sigma[2] > 0 r = r/sigma[2] end;
            end;
        else
            r = x
        end;
        r = min(max(r,-1),1)
        d, a = acos(r), LinRange(0,2*pi,npoints)
        val1, val2 = t*sigma[1]*cos.(a.+(d/2)) .+ centre[1], t*sigma[2]*cos.(a.+(d/2)) .+ centre[2]
        return hcat(val1,val2)
    end;

    function cond_mean_and_cov(X,g;centre=nothing,bary::Bool=false)
        # Filter
        dat = Matrix(filter(x->x[names(X)[1]] == g,X)[:,2:end])

        if centre === nothing
            center = mean(dat,dims=1)
        else
            if size(dat)[2] !== length(centre)
                throw(ArgumentError("Length of centre incorrect"))
            end;
            center = centre[g]
        end;
        
        if size(dat)[1] > 1 
            cov_mat = cov(dat,corrected=true) 
        else
            cov_mat = zeros((2,2))
        end;
        
        if bary cov_mat = cov_mat ./ size(dat)[1] end;
        return center, cov_mat
    end;

    u = sort(unique(X[!,1]))
    stats = Dict(g => cond_mean_and_cov(X,g;centre=centre,bary=bary) for g in u)
    res = reduce(vcat,[ellipse(stats[g][2],centre=stats[g][1],level=level_conf,npoints=npoints) for g in u])
    # Convert to dataframe
    return hcat(DataFrame(names(X)[1] => reduce(vcat,fill.(u,npoints))),DataFrame(res,names(X)[2:end]))
end;
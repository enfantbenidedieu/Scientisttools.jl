# Load packages
using DataFrames
using Statistics
using StatsBase

function conditional_weighted_average(X::AbstractDataFrame,Y::AbstractDataFrame;w=nothing)
    """Compute the conditional weighted average.
    
    """
    n_rows,quant_names,vquant = nrow(X),names(X),Matrix{Float64}(X)
    if w === nothing w = ones(n_rows)/n_rows end;

    function cond_mean(lab)
        vsqual = Y[!,lab]
        modalite = sort(unique(vsqual))
        function mod_mean(mod)
            idx = [i for (i,cat) in enumerate(vsqual) if cat==mod]
            mu = mean(vquant[idx,:],weights(w[idx]),dims=1)
            return hcat(DataFrame("Modalities" => mod),DataFrame(mu,quant_names))
        end;
        return reduce(vcat,[mod_mean(mod) for mod in modalite])
    end;
    return reduce(vcat,[cond_mean(lab) for lab in names(Y)])
end;
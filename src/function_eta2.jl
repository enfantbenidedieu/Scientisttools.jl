# Correlation ratio
using DataFrames

include("get_dummies.jl")

function function_eta2(X::AbstractDataFrame,lab::AbstractString,Y::AbstractDataFrame;w=nothing)

    n_rows,vquant_names,vquant = size(Y)[1],names(Y),Matrix{Float64}(Y)

    # Set weights if noting
    if w === nothing w = ones(n_rows)/n_rows end;

    function fct_eta2(idx)
        tt = get_dummies(X[!,[lab]])
        tt_mat = Matrix{Int64}(tt)
        ni =  sum(tt_mat .* w,dims=1)
        num = sum((tt_mat .* vquant[:,idx]) .* w,dims=1) .^2
        num = sum(num ./ni)
        denom = sum((vquant[:,idx] .^2) .*w)
        return num/denom
    end;
    stats = reduce(hcat,[fct_eta2(i) for i in 1:length(vquant_names)])
    
    if stats isa AbstractFloat
        stats = DataFrame(vquant_names[1] => stats)
    else
        stats = DataFrame(stats,vquant_names)
    end;
    return hcat(DataFrame("Variables" => lab),stats)
end;
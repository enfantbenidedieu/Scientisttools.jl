# Load packages
using StatsBase

function descrstatsw(X,w)
    X = Matrix{Float64}(X)
    n,p = size(X)
    mu = transpose([mean(X[:,j],w) for j in 1:p])
    sigma = [sqrt(sum(w.*(X[:,j].-mean(X[:,j],w)).^2)/sum(w)) for j in 1:p]
    sigma = transpose(sigma)
    return (mu,sigma)
end;
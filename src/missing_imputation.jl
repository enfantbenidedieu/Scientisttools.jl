function missing_imputation(X;method::Symbol=:mean)
    # Make a copy and convert to matrice
    Y = copy(Matrix(X))
    function imputation(x;method::Symbol=:mean)
        if any(i -> isa(i,String), x)
            return coalesce.(x,mode(skipmissing(x))) 
        else
            if method == :mean
                return coalesce.(x,mean(skipmissing(x)))
            elseif  method == :median
                return coalesce.(x,median(skipmissing(x)))
            end;
        end;
    end;

    Y = reduce(hcat,[imputation(Y[:,i],method=method) for i in 1:ncol(X)])

    if X isa AbstractDataFrame
        return DataFrame(Y,names(X))
    else
        return Y
    end;
end;

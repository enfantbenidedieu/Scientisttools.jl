# Load packages
using DataFrames

function reconst(self::NamedTuple;n_components=nothing)

    # Check if self is an object of NamedTuple
    if !(self isa NamedTuple) throw(ArgumentError("self must be NamedTuple object.")) end;

    # Check if PCA or CA model
    if !(self.model in ["pca","ca","mfa"]) throw(ArgumentError("self must be a NamedTuple of PCA, CA or MFA.")) end;

    # Check if n_components
    if n_components !== nothing
        n_components = min(n_components,self.call.n_components)
    else
        n_components = self.call.n_components
    end;

    if self.model == "ca"
        X = Matrix{Int}(self.call.X[!,2:end])
        l, c, total = size(X)[1],size(X)[2],sum(X)
        freq = X ./ total
        col_sum = vec(sum(freq,dims=1))
        row_sum = vec(sum(freq,dims=2))
        eigen_values = vec(self.eig[!,2])
        row_coord = Matrix{Float64}(self.row.coord[!,2:(n_components+1)])
        col_coord = Matrix{Float64}(self.col.coord[!,2:(n_components+1)])
        # Initialisation
        hatX = reshape(row_sum,(l,1)) * reshape(col_sum,(1,c))
        # Add if number of components greater than 0
        if n_components > 0
            U = (row_coord .* sqrt.(row_sum)) ./ transpose(sqrt.(eigen_values[1:n_components]))
            V = (col_coord .* sqrt.(col_sum)) ./ transpose(sqrt.(eigen_values[1:n_components]))
            S = U .* transpose(sqrt.(eigen_values[1:n_components])) * transpose(V)
            hatX = total .* (((S .* sqrt.(row_sum)) .* transpose(sqrt.(col_sum))) .+ hatX)
        end;
        # Convert to data frame
        hatX = hcat(select(self.row.coord, :Modalites),DataFrame(hatX,self.col.coord[!,1]));
    elseif self.model in ["pca","mfa"]
        if self.model == "pca"
            var_coord = Matrix{Float64}(self.var.coord[!,2:(n_components+1)])
        else
            throw(ArgumentError("Not yet implemented"))
        end;
        ind_coord = Matrix{Float64}(self.ind.coord[!,2:(n_components+1)])
        eigen_values = vec(self.eig[!,2])
        hatX = ind_coord * transpose(var_coord ./ transpose(sqrt.(eigen_values[1:n_components])))
        hatX = (hatX .* self.call.std) .+ self.call.means
        # convert to data frame
        hatX = hcat(select(self.ind.coord, :Individuals),DataFrame(hatX,self.var.coord[!,1]));
    else
        throw(ArgumentError("Not yet implemented"))
    end;

    return hatX
end;


using XLSX, DataFrames
#=
include("PCA.jl")
data,columns = XLSX.readtable("./data/decathlon.xlsx","Feuil1");
donnee = DataFrame(data,columns);
res_pca = PCA(donnee;n_components=10);
reconst_pca = reconst(res_pca);
=#
include("CA.jl")
data,columns = XLSX.readtable("./data/children.xlsx","Feuil1");
donnee = DataFrame(data,columns);
res_ca = CA(donnee);
reconst_ca = reconst(res_ca,n_components=1);

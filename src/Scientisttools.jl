module Scientisttools

    include("load_packages.jl")
    include("load_functions.jl")

    export 
    # Principal Components Analysis (PCA) 
    PCA,
    predictPCA,
    supvarPCA,
    get_pca,
    get_pca_ind,
    get_pca_var,
    # Correspondence Analysis (CA)
    CA,
    predictCA,
    supvarCA,
    get_ca,
    get_ca_col,
    get_ca_row,
    # Multiple Correspondence Analysis (MCA)
    MCA,
    predictMCA,
    supvarMCA,
    get_mca,
    get_mca_ind,
    get_mca_var,
    # Factor Analysis of Mixed Data (FAMD)
    FAMD,
    predictFAMD,
    supvarFAMD,
    get_famd,
    get_famd_ind,
    get_famd_var
    # Eigenvalues
    get_eig,
    get_eigenvalue
end 


include("load_packages.jl")
include("load_functions.jl")
include("PCA.jl")
#include("CA.jl")
#include("MCA.jl")
include("FAMD.jl")
include("eigenvalue.jl")
#include("get_pca.jl")
#include("get_ca.jl")
#include("get_mca.jl")
#include("get_famd.jl")
include("datasets.jl")
#include("anova_test.jl")
#include("catdesc.jl")
#include("wcortest.jl")
#include("contdesc.jl")
#include("dimdesc.jl")
include("fviz_famd.jl")


# Principal Components Analysis (PCA)
decathlon = get_dataset("decathlon");
res_pca = PCA(decathlon,ind_sup=42:46,quanti_sup=12:13,quali_sup=14);
#dim_desc_pca = dimdesc(res_pca);
#eig = fviz_screeplot(res_pca,add_labels=true);
#fig = fviz_pca_ind(res_pca);
#habillage = ["course","hauteur","hauteur","hauteur","course","course","hauteur","hauteur","hauteur","course"];
#fig2 = fviz_pca_var(res_pca,habillage=habillage);
#fig3 = fviz_pca_biplot(res_pca);  

# Correspondence Analysis (CA)
#children = get_dataset("children");
#res_ca = CA(children,row_sup=15:18,col_sup=7:9,quali_sup=10);
#fig = fviz_ca_row(res_ca);
#fig2 = fviz_ca_col(res_ca);
#fig3 = fviz_ca_biplot(res_ca);
#dim_desc_ca = dimdesc(res_ca);

# Multiple Correspondence Analysis (MCA)
#races_canines = get_dataset("races_canines");
#res_mca = MCA(races_canines,ind_sup=28:33,quali_sup=8,quanti_sup=9);
#fig = fviz_mca_ind(res_mca,habillage="Fonction");
#fig2 = fviz_mca_mod(res_mca);
#fig3 = fviz_mca_var(res_mca);
#fig4 = fviz_mca_biplot(res_mca); 
#dim_desc_mca = dimdesc(res_mca);

# Factor Analysis of Mixed Data (FAMD)
#autos2005 = get_dataset("autos2005");
#res_famd = FAMD(autos2005,ind_sup=39:45,quanti_sup=14:16,quali_sup=17);
#fig = fviz_famd_ind(res_famd);
#fig2 = fviz_famd_var(res_famd);
#dim_desc_famd = dimdesc(res_famd);


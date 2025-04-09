@testset "PCA" begin
    
    # Create DataFrame
    X = DataFrame((; :Maths => [6.00, 8.00, 6.00, 14.50, 14.00, 11.00, 5.50, 13.00, 9.00],
                     :Physique => [6.00, 8.00, 7.00, 14.50, 14.00, 10.00, 7.00, 12.50, 9.50],
                     Symbol("Français") => [5.00, 8.00, 11.00, 15.50, 12.00, 5.50, 14.00, 8.50, 12.50],
                     :Anglais => [5.50, 8.00, 9.50, 15.00, 12.50, 7.00, 11.50, 9.50, 12.00]))
    
    res = PCA(X,first_col_as_index=false)
    @test res.model == "pca"
    # Test eigenvalue
    @test round.(res.eig[!,2],digits=3) == [2.876, 1.12, 0.004, 0.001]
    # Test individuals factors coordinates
    @test round.(abs.(res.ind.coord[!,2]),digits=3) == [2.743, 1.241, 1.031, 3.138, 2.051, 0.971, 0.335, 0.62, 0.51]
    # Test variables factor coordinates (correlation with axis)
    @test round.(abs.(res.var.coord[!,2]),digits=3) == [0.811, 0.902, 0.753, 0.915]
end
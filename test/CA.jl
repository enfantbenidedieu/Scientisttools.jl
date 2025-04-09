@testset "CA" begin
    
    # Create datafarme.
    X = DataFrame((; :Names => ["Exp.agri", "Patron", "Cadre.sup", "Employé","Ouvrier"],
                     :Droit => [80, 168, 470, 145, 166],
                     :Sciences => [99, 137, 400, 133, 193],
                     :Médecine => [65, 208, 876, 135, 127],
                     :IUT => [58, 62, 79, 54, 129]))

    res = CA(X)
    @test res.model == "ca"
    # Test eigenvalue
    @test round.(res.eig[!,2],digits=3) == [0.082, 0.002, 0.001]
    # Test rows factors coordinates
    @test round.(abs.(res.row.coord[!,2]),digits=3) == [0.410, 0.020, 0.263, 0.142, 0.451]
    # Test rows contributions
    @test round.(abs.(res.row.contrib[!,2]),digits=3) == [16.292, 0.075, 40.401, 3.024, 40.208]
    # Test rows square cosinus
    @test round.(abs.(res.row.cos2[!,2]),digits=3) == [0.987, 0.123, 0.996, 0.670, 0.992]
    # Test columns factor coordinates
    @test round.(abs.(res.col.coord[!,2]),digits=3) == [0.028, 0.16, 0.303, 0.64]
    # Test columns contributions
    @test round.(abs.(res.col.contrib[!,2]),digits=3) == [0.259, 7.945, 41.584, 50.213]
    # Test columns square cosinus
    @test round.(abs.(res.col.cos2[!,2]),digits=3) == [0.165, 0.948, 0.99, 0.989]
end
using Scientisttools
using Test
using DataFrames

tests = ["PCA","CA"]

for test in tests
    include(test*".jl")
end

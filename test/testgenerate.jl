include("testsystem1.jl")
srand(12)
Y, X = sample(20, 100, M)

function testmean(X, μ)
    @test norm(mean(X) - μ) < 1.96*std(X)/sqrt(length(X))
end
# Test the test: [try testmean(randn(1000), 0.0).value catch false end for i in 1:1000]


testmean(X[1, :, 1], 1.0)

function tostatic(M::LinearHomogSystem) 
    d2, d = Kalman.dims(M)
    LinearHomogSystem(SVector{d}(M.x0), 
    SMatrix{d,d}(M.P0),
    SMatrix{d,d}(M.Phi),
    SVector{d}(M.b),
    SMatrix{d,d}(M.Q),
    SVector{d2}(M.y),
    SMatrix{d2,d}(M.H),
    SMatrix{d2,d2}(M.R))
end
M2 = tostatic(M)

include("testsystem.jl")
Y, X = sample(20, 100, M)

# Reproducibility
@test sample(MersenneTwister(1), 5, 5, M) == sample(MersenneTwister(1), 5, 5, M)

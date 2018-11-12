
include("testsystem1.jl")

U = @inferred trace(LinearEvolution(Φ, Gaussian(b, Q)), 0 => Gaussian(x0, P0), endtime(10000))

X = @inferred trace(Sample(LinearEvolution(Φ, Gaussian(b, Q))), 0 => x0, endtime(10000))

@show norm(mean(values(X)) - b) < eps()
@show cov(values(X))
@show cov(values(U)[end])
@test norm(cov(values(X)) - cov(values(U)[end])) < 0.2

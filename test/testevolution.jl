
include("testsystem1.jl")

U = @inferred trace(LinearEvolution(Φ, Gaussian(b, Q)), 0 => Gaussian(x0, P0), endtime(10000))

X = @inferred trace(DynamicIterators.Sample2(LinearEvolution(Φ, Gaussian(b, Q))), 0 => x0, endtime(10000))

@show norm(mean(values(X)) - b) < eps()
@show cov(values(X))
@show cov(values(U)[end])
@test norm(cov(values(X)) - cov(values(U)[end])) < 0.2

using DynamicIterators: dyniterate

Random.seed!(11)
X = @inferred trace(DynamicIterators.Sample2(LinearEvolution(Φ, Gaussian(b, Q))), 0 => x0, endtime(10))
Y = Any[0 => x0]
Random.seed!(11)
Y = trace(LinearEvolution(Φ, Gaussian(b, Q)), Trace(Sample(0 => x0), Y, endtime(10)))
@test X == trajectory(Y)

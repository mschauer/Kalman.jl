include("testsystem.jl")
Random.seed!(11)
Y = Any[]
x = rand(StateObs(Gaussian(x0, P0), M.obs))
u = 1 => x[1]
Y = [1 => x[2]]
Y = trajectory(collecttrace(LinearObservation(LinearEvolution(Φ, Gaussian(b, Q)), H, R), Trace(Sample(u), Y, endtime(10))))
@show Y


Random.seed!(11)
x = rand(StateObs(Gaussian(x0, P0), M.obs))
X = trace(DynamicIterators.Sampled(M), 1 => x, endtime(10))

@test Trajectory(X.t, last.(X.x)) == Y

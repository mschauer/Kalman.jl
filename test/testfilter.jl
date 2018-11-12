
include("testsystem.jl")
Random.seed!(11)
x = rand(StateObs(Gaussian(x0, P0), M.obs))
X = trace(Sample(M), 0 => x, endtime(100))



@show X

Y = collect(t=>y for (t, (x,y)) in pairs(X))
kf = kalmanfilter(M)
@show trace(kf(Y), 0 => Gaussian(x0, P0), endtime(200), register = ((i, x),) -> i > 0)

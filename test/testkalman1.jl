
include("testsystem.jl")
Random.seed!(11)
x = rand(StateObs(Gaussian(x0, P0), M.obs))
X = trace(DynamicIterators.Sampled(M), 1 => x, endtime(100))



@show X

Y = collect(t=>y for (t, (x,y)) in pairs(X))

Xf1, ll1 = kalmanfilter(O, 0 => Gaussian(x0, P0), Y)
@show Xf1

Xf2, ll2 = kalmanfilter(M, 0 => Gaussian(x0, P0), Y)
@show Xf2




kf = kalmanfilter(M)
Xf = trace(kf(Y), 0 => Gaussian(x0, P0), endtime(200), register = ((i, x),) -> i > 0)

@show Xf


@test Xf.t == Xf2.t == Xf1.t
@test first.(Xf.x) == Xf2.x == Xf1.x
@test Xf.x[end][3] == ll2

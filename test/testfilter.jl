
include("testsystem.jl")
Random.seed!(11)
x = rand(StateObs(Gaussian(x0, P0), M.obs))
Y = trace(Sample(M), 0 => x, endtime(100))

@show Y


#kf = KalmanFilter(Y, M)

#@test eltype(kf) === GaussianDistributions.Gaussian{Float64,Float64}
#@test length(kf) == n
#est = trace(kf)

#=
include("testsystem1.jl")
srand(11)
Y, X = sample(20, 100, M)
kalmanfilter(Y, M)

include("testsystem.jl")
srand(11)
Y, X = sample(20, 100, M)
kalmanfilter(Y, M)
=#

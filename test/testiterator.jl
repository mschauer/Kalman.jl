using Base.Test
include("testsystem.jl")
srand(11)
n = 20
Y, X = sample(n, M)
kf = KalmanFilter(Y, M)

@test eltype(kf) === GaussianDistributions.Gaussian{Float64,Float64}
@test length(kf) == n
est = collect(kf)
@test eltype(est) == eltype(kf)

f(t, x, P, Ppred, ll, K) = ll
import Kalman: mappedreturntype
Kalman.mappedreturntype(_, ::Type{typeof(f)}) = Float64
kf2 = Kalman.MappedKalmanFilter(zip(Base.Iterators.countfrom(1),Y), M, f)

@test eltype(kf2) === Float64
@test length(kf2) == n
ll = sum(kf2)
@test typeof(ll) == eltype(kf2)

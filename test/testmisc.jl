using Kalman
using GaussianDistributions
using DynamicIterators

Mo(θ) = let Φ = θ,
    b = 0.0,
    Q = 1.0,

    H = 1.0,
    R = 1.0,

    E = LinearEvolution(Φ, Gaussian(b, Q))

    LinearObservation(E, H, R)
end

θ0 = 1.0
M0 = Mo(θ0)
x = 0.2
Y = [0 => 0.2]

Y = collecttrace(M0.P, Trace(Sample(0 => x), Y, endtime(1000)));


x0 = 1.0
P0 = 1.0

target(x) = DynamicIterators._lastiterate(bind(Y, Mo(x)), Start(Kalman.Filter(0 => Gaussian(x0, P0), 0.0)))
@show target(1.0)
@time target(1.0)
true

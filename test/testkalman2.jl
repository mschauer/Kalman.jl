include("testsystem.jl")
Random.seed!(11)
G0 = Gaussian(x0, P0)
x = rand(StateObs(Φ*G0, M.obs))
X = trace(DynamicIterators.Sampled(M), 1 => x, endtime(3))
@show x

Y = collect(t=>y for (t, (x,y)) in pairs(X))
@show Y

Xf, ll = kalmanfilter(M, 0 => G0, Y)

Xs, ll = rts_smoother(M, 0 => G0, Y)
@show Xs

QL = sqrt(Q)

P0L = sqrt(Φ*P0*Φ' + Q)
x0 = Φ*x0

RL = sqrt(R)

μ = [ x0; Φ*x0; Φ*Φ*x0; x0; Φ*x0; Φ*Φ*x0; ]
L = [    P0L      0    0  0  0  0
       Φ*P0L     QL    0  0  0  0
     Φ*Φ*P0L   Φ*QL   QL  0  0  0
       H*P0L      0    0 RL  0  0
     H*Φ*P0L   H*QL    0  0 RL  0
   H*Φ*Φ*P0L H*Φ*QL H*QL  0  0 RL
]
Σ = L*L'
v = last.(Y)
@show μ
@show Σ
@show v

n = 3
F = [GaussianDistributions.conditional(Gaussian(μ, Σ), 1:i, (1:i).+n, v[1:i]) for i in 1:n]
μf = last.(mean.(F))
Σf = last.(cov.(F))
μs = mean(F[n])
Σs = diag(cov(F[n]))

o = (1:n).+n
ll2 = logpdf(Gaussian(μ[o], Σ[o, o]), v)

@test μf ≈ mean.(values(Xf))
@test Σf ≈ cov.(values(Xf))

@test μs ≈ mean.(values(Xs))
@test Σs ≈ cov.(values(Xs))

@test ll ≈ ll2

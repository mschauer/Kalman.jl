
function backward_rand_kernel(x, Gf, Gpred, Φ) where T
    xf, Pf = meancov(Gf)
    xpred, Ppred = meancov(Gpred)

    J = Pf*Φ'/Ppred
    xs = xf +  J*(x - xpred)
    Ps = Pf - J*Ppred*J'
    rand(Gaussian(xs, Ps))
end

struct BackwardRand{T,S}
    iter::T # iterate backwards through
    M::S
end

function dyniterate(B::BackwardRand, (x, state))
    ϕ = dyniterate(B.iter, state)
    ϕ === nothing && return nothing
    (n, U), state = ϕ
    Gf, Gpred = U[1], U[2]
    Gs = backward_rand_kernel(x, Gf, Gpred, Phi(RTS.M.sys, (n, U)))
    n=>x, (x, state)
end


function backward_sampler(M, prior, Y)
    P = filter(Y, M)

    ϕ = dyniterate(P, Start(prior))
    ϕ === nothing && error("no observations")
    (t, u), state = ϕ

    Xf = trajectory((t => (u[1], u[2]),))
    while true
        ϕ = dyniterate(P, state)
        ϕ === nothing && break
        (t, u), state = ϕ
        push!(Xf, t => (u[1], u[2])) # safe only what is needed
    end

    n = length(Xf)
    U = Xf.x[end]
    Gf, Gpred = U[1], U[2]

    x = rand(Gf)
    xs = fill(NaN*x, n)
    xs[n] = x
    for i in n-1:-1:1
        U = Xf.x[i]
        Gf, Gpredᵒ = U[1], U[2]
        x = backward_rand_kernel(x, Gf, Gpred, Phi(M.sys, i => U))
        xs[i] = x
        Gpred = Gpredᵒ
    end
    Trajectory(Xf.t, xs)
end

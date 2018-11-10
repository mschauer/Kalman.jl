

struct GenericLinearObservation <: AbstractObservation
end
struct GenericLinearEvolution{TPhi,Tb,TQ} <: AbstractEvolution
end



"""
    LinearObservation(H, R)

Observe `y = Hx + v` where ``v ~ N(0, R)``
"""
struct LinearObservation{TH, TR} <: AbstractObservation
    H::TH # dxd
    R::TR # d
end

LinearHomogSystem(G0, Phi, Q, H, R) = LinearStateSpaceModel(G0, LinearEvolution(Phi, Q), LinearObservation(H, R))

function evolve(rng, s, x, t, M::LinearEvolution)
    @assert t-s == 1
    M.Phi*x + rand(rng, M.Q)
end

sample(t, x, M) = sample(Base.GLOBAL_RNG, t, x, M)

function sample(rng::AbstractRNG, t, x, M)
    xx = [x]
    for i in 2:length(t)
        x = evolve(rng, t[i-1], x, t[i], M)
        push!(xx, x)
    end
    t, xx
end

function predict!(s, G::T, t, u, M) where {T}
    if s == t # if `s == t` no time has passed
        G(mean(G) + u, cov(P)), one(M.Phi)
    else
        G(M.Phi*mean(G) + mean(M.Q) + u, M.Phi*cov(G)*M.Phi' + cov(M.Q)), M.Phi
    end
end

function predict!(s, t, U, M::GenericLinearEvolution)
    Phi, Q = U # unpack generic input
    G(Phi*mean(G) + mean(Q), Phi*cov(G)*Phi' + cov(Q)), Phi
end

function observe!(s, t, Y, M::GenericLinearObservation)
    y, H, R = Y  # unpack generic observation
    t, y, H, R
end

function observe!(s, t, Y, MO::LinearObservation)
    t, Y, MO.H, MO.R
end

"""
    kalman_kernel(s, x, P, t, U, Y, SSM) -> t, x, P, Ppred, ll, K

Single Kalman filter step consisting of a prediction step `predict!`, an observation step `observe!`
and a correction step `correct!` with input `U`. Return filtered covariance `P` and predicted `Ppred`

Computes and returns as well the log likelihood of the residual and the Kalman gain.
"""
function kalman_kernel(s, G, t, U, Y, SSM)

    Gpred, Phi, b = predict!(s, G, t, U, SSM)

    t, y, H, R = observe!(s, t, Y, SSM)

    G, yres, S, K = correct!(Gpred, y, H, R, SSM)

    ll = llikelihood(yres, S, SSM)

    t, G, Gpred, ll, K
end

"""
    kalmanfilter!(tt, uu, yy, xx, PP, xxpred, PPpred, M)

Filter obserations in place.
`tt` are `n` (sic!) observation time points, `uu` are `n` control inputs and
`yy` are `n` (generalized) observations (depending on context either vectors or matrices with `size(yy, 2) == n` etc.)
"""
function kalmanfilter!(tt, uu, yy, GG, GGpred, M)
    n = length(tt)

    G = prior(M)

    ll = 0.0
    t = tt[1]

    for i in 1:n
        t, G, Gpred, l, _ = kalman_kernel(t, G, tt[i], uu[.., i], yy[.., i], M)
        GG[i], GGpred[i] = G, Gpred
        ll += l
    end
    GG, GGpred, ll
end

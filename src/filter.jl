

abstract type FilterMethod
end
struct JosephForm <: FilterMethod
end
struct SimpleKalman <: FilterMethod
end

function correct(SSM::LinearStateSpaceModel, u::T, y) where T
    x, Ppred = meancov(u)
    H, R = SSM.obs.H, SSM.obs.R
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S
end

function correct(method::JosephForm, u::T, (v, H)::Tuple{Gaussian}) where T
    x, Ppred = meancov(u)
    y, R = meancov(v)
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S
end

"""

Single Kalman filter with control being the observation, consisting of a prediction step
(evolution of the underlying system) and a correction step `correct`.
Return filtered state `t => (x, P)` and predicted `(xpred, Ppred)`.

Computes and returns as well the log likelihood of the residual.
"""
function dyniterate(M::StateSpaceModel, (s, u)::Pair, (c,)::Control)
    t, y = c
    t, upred = evolve(M.sys, s => u, t)
    u, yres, S = correct(M, upred, y)
    ll = llikelihood(M, yres, S)
    (t => (u, upred, ll)), t => u
end

function dyniterate(M::StateSpaceModel, ::Nothing, (value, c)::NamedTuple{(:value, :control)})
    s, u = value
    t, y = c
    t, upred = evolve(M.sys, s => u, t)
    u, yres, S = correct(M, upred, y)
    ll = llikelihood(M, yres, S)
    (t => (u, upred, ll)), t => u
end

#=
function dyniterate(M::StateSpaceModel, (s, u)::Pair, (c,)::Control)
    t, (v, H) = c
    t, upred = evolve(M.sys, s => u, t)
    u, yres, S, K = correct(JosephForm(), upred, (v, H))
    ll = llikelihood(M, yres, S)
    t => (u, Ppred, ll, K), t => u
end
=#

# fixme: not really clear if being an evolution means that also controlled
# objects are evolutions

function evolve(M::StateSpaceModel, U::Pair, V::Pair)
    _, U = dyniterate(M, U, (control = V,))
    U
end

"""
    kalmanfilter(M, t => x0) -> kf

`kf(iter)` is an iterator over `Gaussian`s or `Distributions` representing
the filtered distribution of `x` where `y` iterates over
(enumerated) signal values.

# Example

```
kf = kalmanfilter(M, 0 => prior) #
est1 = collect(kf(Y1))
est2 = collect(kf(Y2))
```
"""
kalmanfilter(M, prior) = iter -> control(iter, from(M, prior))
kalmanfilter(M) = iter -> control(iter, M)

function kalmanfilter(M, prior, Y)
    P = control(Y, M)

    ϕ = dyniterate(P, nothing, (value=prior,))
    ϕ === nothing && error("no observations")
    (t, u), state = ϕ

    X = trajectory((t => u[1],))
    while true
        ϕ = dyniterate(P, state)
        ϕ === nothing && break
        (t, u), state = ϕ
        push!(X, t => u[1])
    end
    ll = u[3]
    X, ll
end



abstract type FilterMethod
end
struct JosephForm <: FilterMethod
end
struct SimpleKalman <: FilterMethod
end

function correct(method::JosephForm, SSM, u::T, y) where T
    x, Ppred = meancov(u)
    H, R = SSM.obs.H, SSM.obs.R
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S, K
end
function correct(method::JosephForm, u::T, (v, H)::Tuple{Gaussian}) where T
    x, Ppred = meancov(u)
    y, R = meancov(v)
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S, K
end

"""

Single Kalman filter with control being the observation, consisting of a prediction step
(evolution of the underlying system) and a correction step `correct`.
Return filtered state `t => (x, P)` and predicted `(xpred, Ppred)`.

Computes and returns as well the log likelihood of the residual and the Kalman gain.
"""
function dyniterate(M::StateSpaceModel, (s, u)::Pair, (c,)::Control)
    t, y = c
    t, upred = evolve(M.sys, s => u, t)
    u, yres, S, K = correct(JosephForm(), M, upred, y)
    ll = llikelihood(M, yres, S)
    t => (u, upred, ll, K), t => u
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
#=


struct KalmanFilter{T, TM} <: DynamicIterator
    prior::T
    it::Tit
    M::TM
end

#Base.IteratorSize(::KalmanFilter{Tit}) where {Tit} = Base.IteratorSize(Tit) == Base.HasShape() ? Base.HasLength() : Base.iteratorsize(Tit)
#Base.IteratorEltype(::KalmanFilter) = Base.HasEltype()
#Base.eltype(::Type{KalmanFilter{Tit,TM}}) where {Tit,TM} = eltype(TM)

function iterate(kf::KalmanFilter)
    u = kf.prior
    ϕ = iterate(kf.it, st)
    ϕ === nothing && return nothing
    v, state = ϕ
    u = evolve(kf.M, u, v)
    u, (u, state)
end

function iterate(kf::KalmanFilter, state)
    (t, G), st = state
    ϕ = iterate(kf.it, st)
    ϕ === nothing && return nothing
    (tᵒ, y), st2 = ϕ
    u = evolve(kf.M, t => G, tᵒ => y)
    u, (u, st2)
end

length(kf::KalmanFilter) = length(kf.it)
=#

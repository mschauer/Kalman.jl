

abstract type FilterMethod
end
struct JosephForm <: FilterMethod
end
struct SimpleKalman <: FilterMethod
end

function correct!(method::JosephForm, SSM, Gpred::T, y, H, R) where T
    x, Ppred = meancov(Gpred)
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S, K
end

"""
    kalman_kernel(s, G, t, Y, SSM) -> t, G, Ppred, ll, K

Single Kalman filter step consisting of a prediction step `predict!`, an observation step `observe!`
and a correction step `correct!`. Return filtered covariance `P` and predicted `Ppred`

Computes and returns as well the log likelihood of the residual and the Kalman gain.
"""
function kalman_kernel(s, G, t, Y, SSM)

    Gpred, Phi = predict!(s, G, t, SSM)

    t, y, H, R = observe!(s, t, Y, SSM)

    G, yres, S, K = correct!(Gpred, y, H, R, SSM)

    ll = llikelihood(yres, S, SSM)

    t, G, Ppred, ll, K
end


"""
    KalmanFilter(x0, y, M)

Kalman filter as iterator, iterating over `Gaussian`s or `Distributions` representing
the filtered distribution of `x`. Arguments `y` iterates over signal values.

# Example

```
kf = KalmanFilter(Y, M) #
est = collect(kf)
```
"""
struct KalmanFilter{Tit,TM} <: DynamicIterator
    it::Tit
    M::TM
end

#Base.IteratorSize(::KalmanFilter{Tit}) where {Tit} = Base.IteratorSize(Tit) == Base.HasShape() ? Base.HasLength() : Base.iteratorsize(Tit)
#Base.IteratorEltype(::KalmanFilter) = Base.HasEltype()
#Base.eltype(::Type{KalmanFilter{Tit,TM}}) where {Tit,TM} = eltype(TM)


function dyniterate(kf::KalmanFilter)
    G = prior(kf.M)
    ϕ = iterate(kf.it)
    ϕ === nothing && return nothing
    y, st = ϕ
    _, G = kalman_kernel(0.0, G, 0.0, y, kf.M)
    G, (st, G)
end

function dyniterate(kf::KalmanFilter, state)
    st, G = state
    ϕ = iterate(kf.it, st)
    ϕ === nothing && return nothing
    y, st2 = ϕ
    _, G = kalman_kernel(0.0, G, 1.0, y, kf.M)
    G, (st2, G)
end

length(kf::KalmanFilter) = length(kf.it)

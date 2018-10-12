import Base: start, next, done, iteratorsize, iteratoreltype, eltype, length

"""
    KalmanFilter(y, M)

Kalman filter as iterator, iterating over `Gaussian`s representing
the filtered distribution of `x`. Arguments `y` iterates over signal values.

# Example

```
kf = KalmanFilter(Y, M) #
est = collect(kf)
```
"""
struct KalmanFilter{Tit,TM}
    it::Tit
    M::TM
end

Base.iteratorsize(::KalmanFilter{Tit}) where {Tit} = Base.iteratorsize(Tit) == Base.HasShape() ? Base.HasLength() : Base.iteratorsize(Tit) 
Base.iteratoreltype(::KalmanFilter) = Base.HasEltype()
Base.eltype(::Type{KalmanFilter{Tit, T}}) where {Tit, T<:LinearHomogSystem{Tx, TP}} where {Tx, TP} = Gaussian{Tx,TP}

start(kf::KalmanFilter) = start(kf.it), prior(kf.M)
function next(kf::KalmanFilter, state)
    st, N = state
    y, st2 = next(kf.it, st)
    _, x, P = kalman_kernel(0.0, mean(N), cov(N), 1.0, y, kf.M)
    N2 = Gaussian(x, P)
    N2, (st2, N2)
end

done(kf::KalmanFilter, state) = done(kf.it, state[1])

length(kf::KalmanFilter) = length(kf.it)


"""
    MappedKalmanFilter(tuy, M, f)

Kalman filter with model `M` as iterator, iterating over the result `ret` of
`ret = f(t, x, P, Ppred, u, ll, K)`. `tuy` iterates over tuples `(t, u, y)` 
of signal time, control signal and signal value.

Eltype is determined by calling `mappedreturntype(::typeof(f))`
with fallback any.

# Example
```
f(t, x, P, Ppred, u, ll, K) = ll
kf = MappedKalmanFilter(zip(Base.Iterators.countfrom(1), Y), M, f)
ll = sum(kf)
```

For performance, define:
```
import Kalman: mappedreturntype
Kalman.mappedreturntype(_, ::Type{typeof(f)}) = Float64
```
"""
struct MappedKalmanFilter{Tit,TM,F,Tt}
    it::Tit
    M::TM
    f::F
    t0::Tt
end

Base.iteratorsize(::MappedKalmanFilter{Tit}) where {Tit} = Base.iteratorsize(Tit) == Base.HasShape() ? Base.HasLength() : Base.iteratorsize(Tit) 
Base.iteratoreltype(::MappedKalmanFilter) = Base.HasEltype()
Base.eltype(::Type{MappedKalmanFilter{Tit, TM, F}}) where {Tit, TM, F} = mappedreturntype(TM, F)
mappedreturntype(TM, F) = Any

start(kf::MappedKalmanFilter) = start(kf.it), kf.t0, prior(kf.M)
function next(kf::MappedKalmanFilter, state)
    st, s, N = state

    tuy, st2 = next(kf.it, st)
    t, u, y = tuy
    t, x, P, Ppred, ll, K = kalman_kernel(s, mean(N), cov(N), t, u, y, kf.M)  
    ret = kf.f(t, x, P, Ppred, u, ll, K)
    ret, (st2, t, Gaussian(x, P))
end

done(kf::MappedKalmanFilter, state) = done(kf.it, state[1])

length(kf::MappedKalmanFilter) = length(kf.it)
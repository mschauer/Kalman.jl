import Base: iterate, IteratorSize, IteratorEltype, eltype, length

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
struct KalmanFilter{Tit,TM}
    it::Tit
    M::TM
end

Base.IteratorSize(::KalmanFilter{Tit}) where {Tit} = Base.IteratorSize(Tit) == Base.HasShape() ? Base.HasLength() : Base.iteratorsize(Tit) 
Base.IteratorEltype(::KalmanFilter) = Base.HasEltype()
Base.eltype(::Type{KalmanFilter{Tit,TM}}) where {Tit,TM} = eltype(TM)


function iterate(kf::KalmanFilter)
    G = prior(kf.M)
    ϕ = iterate(kf.it)
    ϕ === nothing && return nothing
    y, st = ϕ
    _, G = kalman_kernel(0.0, G, 0.0, y, kf.M)
    G, (st, G)
end

function iterate(kf::KalmanFilter, state)
    st, G = state
    ϕ = iterate(kf.it, st)
    ϕ === nothing && return nothing
    y, st2 = ϕ
    _, G = kalman_kernel(0.0, G, 1.0, y, kf.M)
    G, (st2, G)
end

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

start(kf::MappedKalmanFilter) = 

function iterate(kf::MappedKalmanFilter)
    s, G = kf.t0, prior(kf.M)
    ϕ = iterate(kf.it)
    ϕ === nothing && return nothing
    tuy, st2 = ϕ
    t, u, y = tuy
    t, G, Gpred, ll, K = kalman_kernel(s, G, t, u, y, kf.M)  
    ret = kf.f(t, G, Gpred, u, ll, K)
    ret, (st2, t, G)
end
function next(kf::MappedKalmanFilter, state)
    st, s, G = state
    ϕ = iterate(kf.it, st)
    ϕ === nothing && return nothing
    tuy, st2 = ϕ
    t, u, y = tuy
    t, G, Gpred, ll, K = kalman_kernel(s, G, t, u, y, kf.M)  
    ret = kf.f(t, G, Gpred, u, ll, K)
    ret, (st2, t, G)
end

done(kf::MappedKalmanFilter, state) = done(kf.it, state[1])

length(kf::MappedKalmanFilter) = length(kf.it)
"""
```
LinearStateSpaceModel <: StateSpaceModel

LinearStateSpaceModel(sys, obs, prior)
```

Combines a linear system `sys`, an observations model `obs` and
a `prior` to a linear statespace model in a modular way. See [LinearHomogSystem`](@ref)
for a "batteries included" complete linear system.
"""
struct LinearStateSpaceModel{T,Tsys,Tobs} <: StateSpaceModel
    prior::T
    sys::Tsys
    obs::Tobs
end

prior(M) = M.prior
gausstype(::LinearStateSpaceModel{T}) where {T} = T

predict!(s, G, t, U, M::LinearStateSpaceModel) = predict!(s, G, t, U, M.sys)
predict!(s, G, t, M::LinearStateSpaceModel) = predict!(s, G, t, M.sys)

observe!(s, t, Y, M::LinearStateSpaceModel) = observe!(s, t, Y, M.obs)
prior(M::LinearStateSpaceModel) = M.prior

dims(SSM) = size(SSM.H)

llikelihood(yres, S, SSM) = logpdf(Gaussian(zero(yres), S), yres)

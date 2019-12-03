
abstract type StateSpaceModel <: DynamicIterator
end


struct StateObs{T,S}
    x::T
    obs::S
end

function rand(rng::AbstractRNG, U::StateObs)
    x = rand(rng, U.x)
    y = rand(rng, U.obs(Gaussian(x, false*I)))
    x, y
end

"""
```
LinearStateSpaceModel <: StateSpaceModel

LinearStateSpaceModel(sys, obs)
```

Combines a linear system `sys` and an observations model `obs` and
to a linear statespace model in a modular way.

Evolves `StateObs` objects.
"""
struct LinearStateSpaceModel{Tsys,Tobs} <: StateSpaceModel
    sys::Tsys
    obs::Tobs
end

function evolve(SSM::LinearStateSpaceModel, (t, x)::Pair{<:Any, <:Gaussian})
    t, x = evolve(SSM.sys, t => x)
    t => StateObs(x, SSM.obs)
end

function evolve(SSM::LinearStateSpaceModel, (t, (x,y))::Pair{<:Any, <:Tuple})
    t, x = evolve(SSM.sys, t => x)
    t => StateObs(x, SSM.obs)
end

dims(SSM) = size(SSM.H)

llikelihood(::LinearStateSpaceModel, yres, S) = logpdf(Gaussian(zero(yres), symmetrize!(S)), yres)



abstract type Observation <: DynamicIterator
end

struct GenericLinearObservationModel <: Observation
end

"""
    LinearObservationModel(H, R)

Observe `y = Hx + v` where ``v \\sim N(0, R)``.
"""
struct LinearObservationModel{TH, TR} <: Observation
    H::TH # dxd
    R::TR # d
end

(O::LinearObservationModel)((x, P)::G) where {G} = G(O.H*x, O.H*P*O.H' + O.R)
(O::LinearObservationModel)((x, P)::Gaussian) = Gaussian(O.H*x, O.H*P*O.H' + O.R)

"""
    LinearObservation(P, H, R)

Observe the LinearEvolution `P` using `y = Hx + v`
where ``v \\sim N(0, R)``.

# Examples
```
    O = LinearObservation(LinearEvolution(Φ, Gaussian(b, Q)), H, R)
```
"""
struct LinearObservation{T, TH, TR} <: Observation
    P::T
    H::TH # dxd
    R::TR # d
end


function dyniterate(O::LinearObservation, u)
    (t, x, P), u = @returnnothing dyniterate(O.P, u)
    t => Gaussian(O.H*x, O.H*P*O.H' + O.R), u
end

function dyniterate(O::LinearObservation, u::Sample)
    (t, x), u = @returnnothing dyniterate(O.P, u)
    t => rand(u.rng, Gaussian(O.H*x, O.R)), u
end

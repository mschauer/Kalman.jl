

abstract type Observation <: DynamicIterator
end

struct GenericLinearObservation <: Observation
end

"""
    LinearObservation(H, R)

Observe `y = Hx + v` where ``v ~ N(0, R)``
"""
struct LinearObservation{TH, TR} <: Observation
    H::TH # dxd
    R::TR # d
end

(O::LinearObservation)((x, P)::G) where {G} = G(O.H*x, O.H*P*O.H' + O.R)
(O::LinearObservation)((x, P)::Gaussian) = Gaussian(O.H*x, O.H*P*O.H' + O.R)



function dyniterate(O::LinearObservation, ::Nothing, ((t, (x, P)),)::@NT(control))
    t => Gaussian(O.H*x, O.H*P*O.H' + O.R), nothing
end

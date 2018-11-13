

abstract type Observation <: DynamicIterator
end

struct GenericLinearObservation2 <: Observation
end


struct LinearObservation2{TH, TR} <: Observation
    H::TH # dxd
    R::TR # d
end

(O::LinearObservation2)((x, P)::G) where {G} = G(O.H*x, O.H*P*O.H' + O.R)
(O::LinearObservation2)((x, P)::Gaussian) = Gaussian(O.H*x, O.H*P*O.H' + O.R)

"""
    LinearObservation2(H, R)

Observe `y = Hx + v` where ``v ~ N(0, R)``
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

function dyniterate(E::LinearEvolution, (u, rng)::Sample)
    ϕ = evolve(E, u)
    (t, u) = ϕ
    x = rand(rng, u)
    t => x, Sample(t => x, rng)
end

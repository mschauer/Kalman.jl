module Kalman

using Distributions, GaussianDistributions
using Trajectories, DynamicIterators
using Random, LinearAlgebra

import DynamicIterators: evolve, dyniterate, @returnnothing
import Base: iterate, IteratorSize, IteratorEltype, eltype, length

import Random.rand
import GaussianDistributions: ⊕

meancov(G) = mean(G), cov(G)
meancov(G::Tuple) = G

struct Condition{T} <: Message2
    u::T
    ll::Float64
end
struct Filter{T} <: Message2
    u::T
    ll::Float64
end

macro NT(args...)
    :(NamedTuple{($(args)...,)})
end

include("dyniterate.jl")

export LinearEvolution, Evolution
include("evolution.jl")

export LinearObservation, LinearObservationModel, GenericLinearObservationModel, Observation, Observe
include("observation.jl")

export LinearStateSpaceModel, StateObs
include("statespacemodel.jl")

export kalmanfilter
include("filter.jl")

export rts_smoother
include("smoother.jl")
include("backwardsampler.jl")

include("combinator.jl")

end # module

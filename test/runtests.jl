using Kalman
using Test
using StaticArrays
using Random, LinearAlgebra
using Distributions
using GaussianDistributions
using DynamicIterators
using Trajectories

#include(joinpath("..", "docs", "make.jl"))

Random.seed!(1)

include("testevolution.jl")
include("testobservation.jl")


#include("testgenerate.jl")
include("testfilter.jl")
include("testsmoother.jl")
include("testkalman2.jl")

#include("testiterator.jl")
#include("testtrack.jl")
#include("testinput.jl")

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


include("testkalman1.jl")
include("testkalman2.jl")
include("testsmoother.jl")

include("../example/readme.jl")

# include("testmisc.jl") fixed on master

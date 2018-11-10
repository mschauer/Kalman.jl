using Kalman
using Test
using StaticArrays
using Random, LinearAlgebra
using Distributions
using GaussianDistributions
#include(joinpath("..", "docs", "make.jl"))

# write your own tests here
include("testgenerate.jl")
include("testfilter.jl")
include("testsmoother.jl")
include("testiterator.jl")
include("testtrack.jl")
include("testinput.jl")
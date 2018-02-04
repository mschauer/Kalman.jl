using Kalman
using Base.Test
using StaticArrays

include(joinpath("..", "docs", "make.jl"))


# write your own tests here
include("testgenerate.jl")
include("testfilter.jl")
include("testsmoother.jl")
include("testiterator.jl")
include("testtrack.jl")

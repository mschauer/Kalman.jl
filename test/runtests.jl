using Kalman
using Base.Test
using StaticArrays

# write your own tests here
include("testgenerate.jl")
include("testfilter.jl")
include("testsmoother.jl")

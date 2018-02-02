using Kalman
using Base.Test
using StaticArrays
using LaTeXStrings
include("testsystem1.jl")
srand(11)
using Plots

Y, X = sample(20, 1, M)
Xhat = kalmanrts(Y, M)[2]
scatter(Y[:], label=L"Y (= X_1 = obs)")
plot!(reshape(X,2,20)', color = [:red :blue], label=[L"X_1" L"X_2"] )# label=[L"X_1", L"X_2" ])
plot!(reshape(Xhat,2,20)', line= :dot, color = [:red :blue], label=[L"X_1 (est)" L"X_2 (est)"] )# label=[L"X_1", L"X_2" ])


include("testsystem1.jl")
srand(11)
Y, X = sample(20, 100, M)
kalmanrts(Y, M) 

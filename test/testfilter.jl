include("testsystem1.jl")
srand(11)
Y, X = sample(20, 100, M)
kalmanfilter(Y, M) 

include("testsystem.jl")
srand(11)
Y, X = sample(20, 100, M)
kalmanfilter(Y, M) 
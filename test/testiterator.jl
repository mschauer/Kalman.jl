include("testsystem.jl")
srand(11)
Y, X = sample(20, M)
kf = KalmanFilter(Y, M)

est = collect(kf)
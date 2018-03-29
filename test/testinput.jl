using GaussianDistributions
using Kalman, Base.Test
x0 = 1.0
P0 = 1.0

Phi = 1.0
b = 0.2
Q = 1.0

yshadow = NaN
H = 1.0
R = 1.0
M = LinearHomogSystem(x0, P0, Phi, b, Q, yshadow, H, R)

n = 100

srand(11)
Y, X = sample(n, M)

xxf, PP, PPpred, ll = kalmanfilter!(0:n, Y, copy(X), zeros(n), zeros(n), M) 

prior = Gaussian(Phi*x0, Phi*P0*Phi' + Q)
sys = LinearEvolution(Phi, 0.0, Q)
obs = GenericLinearObservation()
obs2 = LinearObservation(H, R)

M2 = LinearStateSpaceModel(sys, obs, prior)

uu = [b for i in 1:n]

#xxf2, PP2, PPpred2, ll2 = kalmanfilter!(1:n, uu, Y, copy(xxf), copy(PP), copy(PPpred), M2)
xxf2, PP2, PPpred2, ll2 = kalmanfilter!(1:n, uu, [(y, H, R) for y in Y], copy(xxf), copy(PP), copy(PPpred), M2)

@test xxf2 == xxf
@test PP2 == PP
@test PPpred2 == PPpred
@test ll2 == ll


M2 = LinearStateSpaceModel(sys, obs2, prior)


xxf2, PP2, PPpred2, ll2 = kalmanfilter!(1:n, uu, Y, copy(xxf), copy(PP), copy(PPpred), M2)

@test xxf2 == xxf
@test PP2 == PP
@test PPpred2 == PPpred
@test ll2 == ll
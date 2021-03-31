using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕ # independent sum of Gaussian r.v.
using Statistics
using GLMakie
using StaticArrays

function onestep!(p, v)
    v .+= randn(2)
    w = randn(2)
    p .+= v + w
end

x0 = zero(SVector{4, Float64})  #start point
P = SMatrix{4,4,Float64}(I)
# H: observation matrix
H = SMatrix{2,4,Float64}(I) # as only location is observed, not speed
# F: state-trasition
F = SMatrix{4,4,Float64}([
    1 0 1 0
    0 1 0 1
    0 0 1 0
    0 0 0 1
   ])
# Q: the covariance of the process noise
Q = SMatrix{4,4}(Diagonal([0.01, 0.01, 0.1, 0.1]))
# R: the covariance of the observation noise
R = SMatrix{2,2,Float64}(I)


p = zeros(2) # initial poistion
v = ones(2) # initial velocity
n = 100

pp = Gaussian(x0, P)
ps = [pp] # vector of filtered Gaussians
obs = Point2f0[]
for i in 1:n
    global pp
    # predict
    pp = F*pp ⊕ Gaussian(zero(x0), Q) #same as Gaussian(Φ*p.μ, Φ*p.Σ*Φ' + Q)
    # observe
    p, v = onestep(p, v)
    xobs = copy(p)
    push!(obs, xobs)
    # correct
    pp, yres, _ = Kalman.correct(Kalman.JosephForm(), pp, (Gaussian(xobs, R), H))
    push!(ps, pp) # save filtered density
end

path = [Point2f0(mean(p)[1:2]) for p in ps]
vel = [Point2f0(mean(p)[3:4]) for p in ps]

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect())
scatter!(obs)
lines!(obs)
sl = Slider(fig[2, 1], range = 1:n, startvalue = 1)
p = @lift([Point2f0(path[$(sl.value)])])
d = @lift([Point2f0(vel[$(sl.value)])])
arrows!(p, d, arrowcolor = :red, linecolor = :red)#, linewidth = 3, )
w = 0.01only(diff([extrema(vcat(obs...))...]))
on(p) do pp
    pp = p[]
    ps = [pp; pp .+ d[]]
    xlim = extrema(first, ps) .+ (-w, w)
    ylim = extrema(last, ps) .+ (-w, w)
    limits!(ax, xlim..., ylim...)
end
fig




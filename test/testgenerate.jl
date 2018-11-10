include("testsystem1.jl")
Random.seed!(12)

function Kalman.sample(n, t, M::LinearStateSpaceModel)
    x0 = rand(M.prior)
    _, x = sample(t, x0, M.sys)
    y = [H*xᵢ + rand(M.obs.R) for xᵢ in x]
    X = [x]
    Y = [y]
    for i in 1:n
        _, x = sample(t, x0, M.sys)
        y = [H*xᵢ + rand(M.obs.R) for xᵢ in x]
        push!(X, x)
        push!(Y, y)
    end
    Y, X
end

Y, X = sample(20, 1:100, M)

function testmean(X, μ)
    @test norm(mean(X) - μ) < 1.96*std(X)/sqrt(length(X))
end
 
   #= @assert body.head == :block
    stms = Any[]
    stms1 = Any[]#

 
    for stm in body.args
        dump(stm)
        if stm isa Expr
            if stm.head == :if && stm.args[1] isa Expr && stm.args[1].head == :$ && stm.args[1].args[1] == :first
                push!(stms1, stm.args[2])
                push!(stms, stm.args[3])
                println("hea")
            else
                push!(stms1, stm)
                push!(stms, stm)
            end
        else
            push!(stms1, stm)
            push!(stms, stm)
        end
    end
    =#


tostatic(::Val{d}, G::Gaussian) where {d} = Gaussian(SVector{d}(mean(G)), SMatrix{d,d}(cov(G)))
tostatic(::Val{d}, G::Matrix) where {d} = SMatrix{d,d}(G)

function tostatic(M) 
    d2, d = Kalman.dims(M)
    LinearHomogSystem(tostatic(Val(d), M.prior),
    SMatrix{d,d}(M.sys.Phi),
    SVector{d}(M.sys.b),
    SMatrix{d,d}(M.sys.Q),
    SMatrix{d2,d}(M.obs.H),
    tostatic(Val(d), M.obs.R))
end
M2 = tostatic(M)

include("testsystem.jl")
Y, X = sample(20, 100, M)

# Reproducibility
@test sample(MersenneTwister(1), 5, 5, M) == sample(MersenneTwister(1), 5, 5, M)

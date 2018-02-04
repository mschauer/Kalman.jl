"""
    track(A::Matrix{Float64}, p, h = 10) -> p2, err

Track blurry lightsource by applying a window with half-width `h` at 
an approximate location `p = (i,j)` and find the
average weighted location of points with *high* light intensity. 

Gives an standard error estimate.
"""
function track(img::Matrix, p, h = 10)
    i, j = round.(Int, p)
    CR = CartesianRange(CartesianIndex(i - h, j - h), CartesianIndex(i + h, j + h))
    μ = mean(img[ci] for ci in CR)
    C = sum(max(img[ci] - μ, 0) for ci in CR) 

    xhat = sum(max(img[ci] - μ, 0)*ci[1] for ci in CR)/C
    yhat = sum(max(img[ci] - μ, 0)*ci[2] for ci in CR)/C

    xerr =  sum(max(img[ci] - μ, 0)*(ci[1] - xhat)^2 for ci in CR)/C
    yerr =  sum(max(img[ci] - μ, 0)*(ci[2] - yhat)^2 for ci in CR)/C

    err = sqrt(xerr + yerr)
    (xhat, yhat), err
end
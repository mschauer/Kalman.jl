A = zeros(10,10)
A[3,4] = 1
A[5,6] = 1
A[10,1] = 10
p, err = track(A, (5,5), 3)
@test norm(p .- [4,5]) < sqrt(eps())
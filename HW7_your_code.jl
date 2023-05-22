
#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix A and returns the 
# kmax × kmax supmatrix of itsupper Hessenberg form 
# The matrix A should be only accessed through 
# (kmax - 1) matrix-vector products
function arnoldi(A, q1, kmax)
    #From Darve and Wooters
        m = size(A)[1]
        H = zeros(kmax, kmax)
        q = zeros(m, kmax)
        q[:,1] = q1 / (q1' * q1)^0.5
        
        # Arnoldi iteration
        for k = 1:kmax
            if k > 1 
                q[:,k] = q1 / H[k,k-1]
            end
            q1 = A * q[:,k] 
            for i=1:k 
                H[i,k] = q[:,i]' * q1
                q1 -= H[i,k] * q[:,i]
            end
            if k<kmax
                H[k+1,k] = (q1' * q1)^0.5
            end
        end

        return H
end
    


#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and returns the 
# kmax × kmax supmatrix of its tridiagonal form 
# computed by the Lanczos iteration.
# The matrix A should be only accessed through 
# (kmax - 1) matrix-vector products
# The output vectors should be the diagonal (α)
# and the offdiagonal (β) of the tridiagonal matrix
function lanczos(A, q1, kmax)
    # From Darve and Wooters
    m = size(A)[1]
    T = zeros(kmax, kmax)
    Q = zeros(m, kmax)
    beta = norm(q1)
    q0 = zeros(m)

    # Lanczos iteration
    for k = 1:kmax
        Q[:, k] = q1 / beta
        q1 = A * Q[:, k]
        alpha = Q[:, k]'*q1
        T[k, k] = alpha
        if k > 1
            T[k - 1, k] = beta
            T[k, k - 1] = beta
        end
        q1 = q1 - alpha * Q[:, k] - beta * q0
        q0 = Q[:, k]
        beta = norm(q1)
    end

    α = diag(T)
    β = diag(T, 1)
    return α, β
end

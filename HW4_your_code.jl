#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix A and returns 
# a reduced QR factorization with factors Q and R.
# It should not modify A
function classical_gram_schmidt(A)
    # YOUR CODE HERE
    # Algorithm 7.1 from Trefethen and Bau
    m, n = size(A)
    Q = zeros(m, n)
    R = zeros(n, n)
    for j = 1:n
        v = A[:,j]
        for i = 1:j-1
            R[i,j] = Q[:,i]'*A[:,j]
            v -= R[i,j]*Q[:,i]
        end
        R[j,j] = norm(v)
        Q[:,j] = v/R[j,j]
    end
    return Q, R
end

#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and returns 
# a reduced QR factorization with factors Q and R.
# It should not modify A
function modified_gram_schmidt(A)
    # YOUR CODE HERE
    # Algorithm 8.1 from Trefethen and Bau
    m, n = size(A)
    Q = zeros(m, n)
    R = zeros(n, n)
    V = copy(A)
    for j = 1:n
        R[j,j] = norm(V[:,j])
        Q[:,j] = V[:,j]/R[j,j]
        for i = j+1:n
            R[j,i] = Q[:,j]'*V[:,i]
            V[:,i] -= R[j,i]*Q[:,j]
        end
    end
    return Q, R
end

#----------------------------------------
# Problem c
#----------------------------------------
# This function takes in a matrix A 
# and computes its QR factorization in place,
# using householder reflections.
# It should not allocate any memory. 
@views function householder_QR!(A)
    # YOUR CODE HERE
    m, n = size(A)
    for k = 1:n 
        mod_x = sign(A[k:m, k][1])*(A[k:m, k]'*A[k:m, k])^0.5
        A[k:m, k][1] += mod_x
        v1 = A[k:m, k][1]
        for i = k:m
            A[i, k] /= v1
        end
        for i = k:n
            vnorm = (A[k:m, k]' * A[k:m, k])
            vA = A[k:m, k]' * A[k:m, i]
            for j = k:m
                A[j, i] -= 2 * A[j, k] * vA/vnorm
            end
        end
        A[k,k] = -1*mod_x
    end
    return A  
end

#----------------------------------------
# Problem d
#----------------------------------------
# These two functions take in the housholder
# QR factorization from part c and multiply them
# to a vector (mul) or solve the least squares 
# problem in A (div), in place.
# They should not allocate any memory and instead
# use the preallocated output vector to record the result. 
@views function householder_QR_mul!(out_mul, x, QR)
    # YOUR CODE HERE
    m, n = size(QR)
    for i = 1:length(out_mul)
        out_mul[i] = 0
    end
    # Calculating Rx
    for i = 1:m
        for j = i:n
            out_mul[i] += QR[i,j]*x[j]
        end
    end
    # Calculation QRx
    for k = n:-1:1
        v = QR[k:m, k]
        for i = 1:length(v)
            v[i] *= -1
        end
        v[1] = 1
        vnorm = v' * v
        vscal = v' * out_mul[k:m]
        for i = k:m
            out_mul[i] -= 2 * QR[i,k] * vscal/vnorm
        end
    end
end

function householder_QR_div!(out_div, b, QR)
    # YOUR CODE HERE
    # Calculating Q*b
    m, n = size(QR)
    for k = 1:n
        v = QR[k:m, k]
        for i = 1:length(v)
            v[i] *= -1
        end
        v[1] = 1
        vnorm = v' * v
        vscal = v' * b[k:m]
        for i = k:m
            b[i] -= 2 * QR[i,k] * vscal/vnorm
        end
    end
    out_div = b[1:n]
    U = triu(QR)[1:n,:]
    # back substitution to solve Ux = b
    for i = n:-1:1
        for j = i+1:n
            out_div[i] -= U[i,j] * out_div[j]
        end
        out_div[i] /= U[i,i]
    end
end
#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in the LU factorization
# together with the permutation P (in the form of
# an array of integers such that P[i] = j means that multiplication with P moves the i=th row to the j-th position)
# It should modify the input variable x in place
function substitution!(x, LU, P = 1 : size(x, 1)) 
   # YOUR CODE HERE
    permute!(x, P)
    # forward substitution
    for i = 2 : size(x)[1]
        for j = 1 : i-1
            x[i] = x[i] - LU[i,j] * x[j]
        end
    end
    # backward substitution
    for j = size(x)[1] : -1 : 1
        for i = j+1 : size(x)[1]
            x[j] = x[j] - LU[j,i] * x[i]
        end
        x[j] = x[j] / LU[j,j]
    end
end

#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and modifies
# it in place, such that is contains the stricly
# lower triangular part of L together with the
# upper triangular part of U.
function unpivoted_LU!(A)
    # YOUR CODE HERE
    #Algorithm 20.1 from Trefethen and Bau
    m = size(A)[1]
    # Loop over the columns
    for j = 1:m
        # Loop over the rows
        for i = j+1:m
            A[i,j] = A[i,j]/A[j,j]
            for k = j+1:m
                A[i,k] = A[i,k] - A[i,j]*A[j,k]
            end
        end
    end
end

#----------------------------------------
# Problem d
#----------------------------------------
# This function takes in a matrix A and modifies
# it in place, such that is contains the stricly
# lower triangular part of L together with the
# upper triangular part of U.
# It uses row-pivoting and stores the resulting 
# row permutation in the array P
function pivoted_LU!(A)
    # The array that will be used to keep track of the permutation
    P = collect(1 : size(A, 1))
    # YOUR CODE HERE
    # Algorithm 21.1 from Trefethen and Bau
    n = size(A)[1]
    for k = 1:n-1
        # Finding the pivot row p for column k
        p = argmax(abs.(A[k:n,k])) + k - 1
        if p != k
            # Swap rows k and p in A and update P
            A[k,:], A[p,:] = A[p,:], A[k,:]
            P[k], P[p] = P[p], P[k]
        end
        # Finding multipliers
        for i = k+1:n
            A[i,k] = A[i,k] / A[k,k]
        end
        # Updating submatrix A[k+1:n,k+1:n]
        for i = k+1:n
            for j = k+1:n
                A[i,j] = A[i,j] - A[i,k] * A[k,j]
            end
        end
    end
    return P
end

#----------------------------------------
# Problem e
#----------------------------------------
# Creates an m × m matrix with a particularly 
# large growth factor
function growth_matrix(m)
    # YOUR CODE GOES HERE
    A = zeros(m, m)
    for i = 1:m
        A[i,m] = 1
    end
    for i = 1:m
        for j = 1:m
            if i == j
                A[i,j] = 1
            elseif i>j
                A[i,j] = -1
            end
        end
    end
    return A
end

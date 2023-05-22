# This function takes in a matrix A and a vector v and writes their product into the vector v
function u_is_A_times_v!(u, A, v)
    # Write your function here! 

    for i in 1:size(u)[1]
        u[i] = 0
    end

    for j in 1:size(A)[2]
        for i in 1:size(A)[1]
            u[i] +=  A[i,j] * v[j]
        end
   end

end

# This function takes in matrices ABC and writes B times C into the matrix A
function A_is_B_times_C!(A, B, C)
    # Write your function here! 

    for i in 1:size(A)[1]
        for j in 1:size(A)[2]
            A[i,j] = 0
        end
    end

    for j in 1:size(C)[2]
        for k in 1:size(C)[1]
            for i in 1:size(B)[1]
                A[i,j] += B[i,k]*C[k,j]
            end
        end
    end

end


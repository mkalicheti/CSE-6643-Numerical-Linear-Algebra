# This macro helps optimize the innermost loop in the kernel on your machine. You should not need to use it anywhere else
using LoopVectorization: @turbo
# This function takes in matrices ABC and adds B times C to the matrix A
function add_to_A_B_times_C!(A, B, C)
    @turbo for j in axes(C, 2)
        for k in axes(B, 2)
            for i in axes(A, 1)
                A[i, j] += B[i, k] * C[k, j]
            end
        end
    end
end

# This function takes in matrices ABC and adds B times C to the matrix A
# It uses blocking into blocks of size bks
# Make sure that your function does not allocate memory
function add_to_A_B_times_C!(A, B, C, bks)
    for jj in 1:bks:size(C, 2)
        for kk in 1:bks:size(B, 2)
            for ii in 1:bks:size(A, 1)
                @inbounds add_to_A_B_times_C!(@view(A[ii:min(ii + bks - 1, size(A, 1)), jj:min(jj + bks - 1, size(C, 2))]), @view(B[ii:min(ii + bks - 1, size(A, 1)), kk:min(kk + bks - 1, size(B, 2))]), @view(C[kk:min(kk + bks - 1, size(B, 2)), jj:min(jj + bks - 1, size(C, 2))]))
            end
        end
    end
end

# Implements a recursive, cache oblivious algorithm
# complete this skeleton
function oblivious_add_to_A_B_times_C!(A, B, C, bks)
    i_size = size(A, 1)
    j_size = size(C, 2)
    k_size = size(B, 2)

    # If we want to further subdivide
    if min(i_size, j_size, k_size) > bks
        i_mid = i_size ÷ 2
        j_mid = j_size ÷ 2
        k_mid = k_size ÷ 2
    
        # Recursively solve the 8 subproblems
        oblivious_add_to_A_B_times_C!(@view(A[1:i_mid, 1:j_mid]), @view(B[1:i_mid, 1:k_mid]), @view(C[1:k_mid, 1:j_mid]), bks)
        oblivious_add_to_A_B_times_C!(@view(A[1:i_mid, 1:j_mid]), @view(B[1:i_mid, k_mid+1:k_size]), @view(C[k_mid+1:k_size, 1:j_mid]), bks)
        oblivious_add_to_A_B_times_C!(@view(A[i_mid+1:i_size, 1:j_mid]), @view(B[i_mid+1:i_size, 1:k_mid]), @view(C[1:k_mid, 1:j_mid]), bks)
        oblivious_add_to_A_B_times_C!(@view(A[i_mid+1:i_size, 1:j_mid]), @view(B[i_mid+1:i_size, k_mid+1:k_size]), @view(C[k_mid+1:k_size, 1:j_mid]), bks)
        oblivious_add_to_A_B_times_C!(@view(A[1:i_mid, j_mid+1:j_size]), @view(B[1:i_mid, 1:k_mid]), @view(C[1:k_mid, j_mid+1:j_size]), bks)
        oblivious_add_to_A_B_times_C!(@view(A[1:i_mid, j_mid+1:j_size]), @view(B[1:i_mid, k_mid+1:k_size]), @view(C[k_mid+1:k_size, j_mid+1:j_size]), bks)
        oblivious_add_to_A_B_times_C!(@view(A[i_mid+1:i_size, j_mid+1:j_size]), @view(B[i_mid+1:i_size, 1:k_mid]), @view(C[1:k_mid, j_mid+1:j_size]), bks)
        oblivious_add_to_A_B_times_C!(@view(A[i_mid+1:i_size, j_mid+1:j_size]), @view(B[i_mid+1:i_size, k_mid+1:k_size]), @view(C[k_mid+1:k_size, j_mid+1:j_size]), bks)
    else
        add_to_A_B_times_C!(A, B, C)
    end
end


# #----------------------------------------
# # Problem a
# #----------------------------------------
# # This function takes in a matrix T and modifies it 
# # in place to Hessenberg form using Householder reduction.
@views function hessenberg_form!(A)
    m = size(A, 1)
    #Algorithm 26.1 from Trefethen and Bau
    for k = 1:(m-2)
        x = A[k+1:m, k]
        v_k = sign(x[1]) * norm(x) * [1.0; zeros(m-k-1)] + x
        v_k = v_k / norm(v_k)
        A[(k+1):m, k:m] -= 2 * v_k * (v_k' * A[(k+1):m, k:m])
        A[1:m, (k+1):m] -= 2 * (A[1:m, (k+1):m] * v_k) * v_k'
    end
    #setting the elements below subdiagonal to 0; they are already nearly o
    for i = 3:m
        for j in 1:(i-2)
            A[i,j] = 0
        end
    end
end


#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix T in Hessenberg form
# and runs a single iteration of the unshifted QR Algorithm 
# using Givens rotations

#-----helper functions-----
# This function takes in two numbers a and b and returns
# the cosine and sine of the Givens rotation that zeros out b
function givens_rotation(a, b)
    if b == 0
        c = 1
        s = 0
    else
        if abs(b) > abs(a)
            r = a / b
            s = 1 / sqrt(1 + r^2)
            c = s * r
        else
            r = b / a
            c = 1 / sqrt(1 + r^2)
            s = c * r
        end
    end
    return c, s
end

# This function takes in a matrix A and returns the
# matrix Q and R such that A = QR
# using Givens rotations
function qrgivens!(A)
    m, n = size(A)
    Q = Matrix{Float64}(I, m, m)

    for j = 1:n
        for i = m:-1:(j + 1)
            G = Matrix{Float64}(I, m, m)
            c, s = givens_rotation(A[i - 1, j], A[i, j])
            G[[i - 1, i], [i - 1, i]] = [c -s; s c]
            A .= G' * A
            Q .= Q * G
        end
    end
    return Q, A
end

#--- main function ---
# This function takes in a matrix T in Hessenberg form
# and runs a single iteration of the unshifted QR Algorithm
# using Givens rotations
function givens_qr!(T)
    m, n = size(T)
    Q, R = qrgivens!(T)
    T .= R * Q
end


#----------------------------------------
# Problem c
#----------------------------------------
# This function takes in a matrix T in Hessenberg form and 
# implements the practical QR algorithm with shifts. 
# The input shift dictates which shift type your 
# algorithm should use. For shift = "single" implement the single shift 
# and for shift = "wilkinson" implement the Wilkinson shift

function practical_QR_with_shifts!(T, shift, eigPracQR = [], mu_storage = [])
    m, n = size(T)
    tol = 1e-8
    max_iter = 1000
    for iter = 1:max_iter
        for k = (n-1):-1:1
            if abs(T[k+1, k]) < tol * (abs(T[k, k]) + abs(T[k+1, k+1]))
                T[k+1, k] = 0
            end
        end

        is_converged = true
        for k = 1:(n-1)
            if T[k+1, k] != 0
                is_converged = false
                break
            end
        end

        if is_converged
            break
        end
        # finding the shift
        if shift == "single"
            mu = T[n, n]
        elseif shift == "wilkinson"
            d = (T[n-1, n-1] - T[n, n]) / 2
            sign_d = sign(d) == 0 ? 1 : sign(d)
            mu = T[n, n] - sign_d * T[n, n-1]^2 / (abs(d) + sqrt(d^2 + T[n, n-1]^2))
        else
            error("Invalid shift type")
        end

        # QR iteration on shifted T
        push!(mu_storage, mu) # for 4d
        T .-= mu * Matrix{Float64}(I, m, m)
        givens_qr!(T)
        T .+= mu * Matrix{Float64}(I, m, m)
        push!(eigPracQR, maximum(diag(T))) # for 4d
    end
    return T
end

#----------------------------------------
# Problem d
#----------------------------------------



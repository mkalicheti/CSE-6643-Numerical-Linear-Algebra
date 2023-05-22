import Pkg.instantiate
instantiate()
using LinearAlgebra
using CairoMakie

#-------------------
#-----Functions-----
#-------------------

#---3(a)---
function eigvalmatrix(A, m)
    # Fill A with 0s
    for i in 1:m
        for j in 1:m
            A[i, j] = 0
        end
    end
    # Fill in eigenvalues along the diagonal
    for i in 1:m
        A[i, i] = 3^(m-i)
    end
    # Find a nonsingular matrix P
    P = rand(m, m)
    while det(P) == 0
        P = rand(m, m)
    end
    # Find the matrix A = PAP^-1
    A = P * A * inv(P)
    return A
end

#---3(b), 3(c), 3(d)---
function orthogonal_iteration(A, Q, k_max, diag_step, block2norm)
    AA = copy(A)
    # Initialize Q as the identity matrix
    for i in 1:size(A)[1]
        for j in 1:size(A)[1]
            if i == j
                Q[i, j] = 1
            else
                Q[i, j] = 0
            end
        end
    end
    # Perform the orthogonal iteration
    for k in 1:k_max
        Z = A * Q
        F = qr(Z)
        Q = F.Q
        println("\nDiagonal of R at iteration $k:")
        for j = 1:size(A)[1]
            #---3(b)---
            print(round(F.R[j, j], digits=4), " ")
            #---3(c)---
            diag_step[k, j] = F.R[j, j]
        end
        #---3(d)---
        p = 4
        Ak = Q' * AA * Q
        block2norm[k] = norm(Ak[p+1:end, 1:p])
    end
    return Q
end

#---3(c)---
# Using the 3 equations in the question
function theoretical_convergence(eigenvalues, p, k)
    if p == 1
        return abs(eigenvalues[p+1] / eigenvalues[p])^k
    elseif p == length(eigenvalues)
        return abs(eigenvalues[p] / eigenvalues[p-1])^k
    else
        return max(abs(eigenvalues[p+1] / eigenvalues[p])^k, abs(eigenvalues[p] / eigenvalues[p-1])^k)
    end
end

#-------------------
#----Driver Code----
#-------------------

m = 8
A = randn(8,8)
outa = eigvalmatrix(A, m) # 3(a)

Q = randn(8,8)
k_max = 5 # max iters
diag_step = randn(k_max, m) # 3(c)
block2norm = randn(k_max) # 3(d)

X = orthogonal_iteration(outa, Q, k_max, diag_step, block2norm) # 3(b), 3(c), 3(d)
println()
println()

# 3(c)
# Calculate the true eigenvalues
eigenvalues = [3^(i-1) for i in m:-1:1]

convergence_1 = diag_step[:, 1]
convergence_2 = diag_step[:, 2]
convergence_3 = diag_step[:, 3]

# Normalise the convergence rates
for i in 1:k_max
    convergence_1[i] = abs(convergence_1[i] - eigenvalues[1])/eigenvalues[1]
    convergence_2[i] = abs(convergence_2[i] - eigenvalues[2])/eigenvalues[2]
    convergence_3[i] = abs(convergence_3[i] - eigenvalues[3])/eigenvalues[3]
end

# Calculate the theoretical convergence rates
theoretical_1 = [theoretical_convergence(eigenvalues, 1, k) for k in 1:k_max]
theoretical_2 = [theoretical_convergence(eigenvalues, 2, k) for k in 1:k_max]
theoretical_3 = [theoretical_convergence(eigenvalues, 3, k) for k in 1:k_max]

# Transform the data using log10
log_convergence_1 = log10.(convergence_1)
log_convergence_2 = log10.(convergence_2)
log_convergence_3 = log10.(convergence_3)
log_theoretical_1 = log10.(theoretical_1)
log_theoretical_2 = log10.(theoretical_2)
log_theoretical_3 = log10.(theoretical_3)

#log plots
figure1 = Figure(resolution = (800, 600), font = "sans")
axis1 = Axis(figure1[1, 1], title = "Convergence p=1", ylabel = "log 10", xlabel = "Iteration")
lines!(axis1, 1:k_max, log_convergence_1, color = :blue, label = "log Experimental")
lines!(axis1, 1:k_max, log_theoretical_1, linestyle = :dash, color = :red, label = "Theoretical")
axislegend(axis1)

axis2 = Axis(figure1[2, 1], title = "Convergence p=2", ylabel = "log 10", xlabel = "Iteration")
lines!(axis2, 1:k_max, log_convergence_2, color = :blue, label = "log Experimental")
lines!(axis2, 1:k_max, log_theoretical_2, linestyle = :dash, color = :red, label = "Theoretical")
axislegend(axis2)

axis3 = Axis(figure1[3, 1], title = "Convergence p=3", ylabel = "log 10", xlabel = "Iteration")
lines!(axis3, 1:k_max, log_convergence_3, color = :blue, label = "log Experimental")
lines!(axis3, 1:k_max, log_theoretical_3, linestyle = :dash, color = :red, label = "Theoretical")
axislegend(axis3)

#abs plots
figure2 = Figure(resolution = (800, 600), font = "sans")
axis1 = Axis(figure2[1, 1], title = "Convergence p=1", ylabel = "Absolute", xlabel = "Iteration")
lines!(axis1, 1:k_max, convergence_1, color = :blue, label = "Absolute Experimental")
lines!(axis1, 1:k_max, theoretical_1, linestyle = :dash, color = :red, label = "Theoretical")
axislegend(axis1)

axis2 = Axis(figure2[2, 1], title = "Convergence p=2", ylabel = "Absolute", xlabel = "Iteration")
lines!(axis2, 1:k_max, convergence_2, color = :blue, label = "Absolute Experimental")
lines!(axis2, 1:k_max, theoretical_2, linestyle = :dash, color = :red, label = "Theoretical")
axislegend(axis2)

axis3 = Axis(figure2[3, 1], title = "Convergence p=3", ylabel = "Absolute", xlabel = "Iteration")
lines!(axis3, 1:k_max, convergence_3, color = :blue, label = "Absolute Experimental")
lines!(axis3, 1:k_max, theoretical_3, linestyle = :dash, color = :red, label = "Theoretical")
axislegend(axis3)

save("plot3c_semilog.png", figure1)
save("plot3c_abs.png", figure2)

# 3(d)
analytical_estimate = [abs(eigenvalues[5]/eigenvalues[4])^k for k in 1:k_max]

figure3 = Figure(resolution = (800, 600), font = "sans")
axis1 = Axis(figure3[1, 1], title = "3(d)", xlabel = "Iteration")
line1 = lines!(axis1, 1:k_max, block2norm, color = :blue, label = "2-norm of Ak")
line2 = lines!(axis1, 1:k_max, analytical_estimate, linestyle = :dash, color = :red, label = "Analytical estimate")
axislegend(axis1)
axis2 = Axis(figure3[2, 1], title = "3(d) semilog", xlabel = "Iteration")
line3 = lines!(axis2, 1:k_max, log10.(block2norm), color = :blue, label = "2-norm of Ak")
line4 = lines!(axis2, 1:k_max, log10.(analytical_estimate), linestyle = :dash, color = :red, label = "Analytical estimate")
axislegend(axis2)
save("plot3d.png", figure3)
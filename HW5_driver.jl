import Pkg.instantiate
instantiate()
using BenchmarkTools: @ballocated
using LinearAlgebra
using CairoMakie
include("HW5_your_code.jl")

m = 8
A = randn(8,8)
A_copy = copy(A)
outa = eigvalmatrix(A, m)
outa_copy = copy(outa)  

Q = randn(8,8)
k_max = 5 #max iters
diag_step = randn(k_max, m)
block2norm = randn(k_max)

#3(c)
# Calculate the true eigenvalues
eigenvalues = [3^(i-1) for i in m:-1:1]

X = orthogonal_iteration(outa, Q, k_max, diag_step, block2norm)
convergence_1 = diag_step[:, 1]
convergence_2 = diag_step[:, 2]
convergence_3 = diag_step[:, 3]
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
# plotting
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
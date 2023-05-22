using BenchmarkTools: @ballocated
using LinearAlgebra: I, norm, triu, tril, tr, diagm, diag, qr, eigvals, istriu
using CairoMakie
include("HW7_your_code.jl")


#----------------------------------------
# Problem a 
#----------------------------------------
########################################
m = 5
A = rand(m,m)
A = A * A' 
H = arnoldi(A, randn(m), m)
@assert istriu(H[2 : end, 1 : (end - 1)])
@assert eigvals(H) ≈ eigvals(A) 

println("Passed part (a) test")


#----------------------------------------
# Problem b 
#----------------------------------------
########################################
m = 5
A = rand(m,m)
A = A * A' 
α, β = lanczos(A, randn(m), m)
H = diagm(-1 => β, 0 => α, 1 => β)
@assert eigvals(H) ≈ eigvals(A) 
println("Passed part (b) test")

#----------------------------------------
# Problem c
#----------------------------------------
########################################

# Your code here

# Helper function to generate symmetric matrices
function generate_symmetric_matrix(m)
    A = rand(m, m)
    return A * A'
end

# Error metric
function eigenvalue_error(true_eigvals, approx_eigvals)
    return norm(true_eigvals - approx_eigvals)
end

# Test matrices
matrix_sizes = [5, 10, 20, 50]
kmax_values = 1:50
num_matrices = length(matrix_sizes)

# Preallocate error storage
arnoldi_errors = zeros(num_matrices, length(kmax_values))
lanczos_errors = zeros(num_matrices, length(kmax_values))

# Perform analysis
for (i, m) in enumerate(matrix_sizes)
    A = generate_symmetric_matrix(m)
    true_eigvals = sort(eigvals(A))
    
    for (j, kmax) in enumerate(kmax_values)
        H_arnoldi = arnoldi(A, randn(m), kmax)
        approx_eigvals_arnoldi = sort(eigvals(H_arnoldi), by=real)
        arnoldi_errors[i, j] = eigenvalue_error(true_eigvals[1:min(kmax, m)], approx_eigvals_arnoldi[1:min(kmax, m)])
        
        α, β = lanczos(A, randn(m), kmax)
        H_lanczos = diagm(-1 => β, 0 => α, 1 => β)
        approx_eigvals_lanczos = sort(eigvals(H_lanczos), by=real)
        lanczos_errors[i, j] = eigenvalue_error(true_eigvals[1:min(kmax, m)], approx_eigvals_lanczos[1:min(kmax, m)])
    end
end

# Plot error vs kmax
figure = Figure(resolution=(1000, 800))

for i in 1:num_matrices
    ax = Axis(figure[i, :], xlabel="kmax", ylabel="Error", title="Matrix size = $(matrix_sizes[i])")
    lines!(ax, kmax_values, arnoldi_errors[i, :], label="Arnoldi")
    lines!(ax, kmax_values, lanczos_errors[i, :], label="Lanczos")
    axislegend(ax)
end

figure
save("error_vs_kmax.png", figure)

#Darve and Wooters Exercise 6.3
n = 128
D = cos.(LinRange(0,π,n))
D[div(n,2)] = 2
Q, = qr(rand(n,n))
x = LinRange(0,1,n)
A_ghost = Q * diagm(0 => D) * Q'
q = zeros(n,n)
true_eigvals_ghost = sort(eigvals(A_ghost))

# Random starting vector
q[:,1] = rand(n)
q[:,1] /= norm(q[:,1])

α, β = lanczos(A_ghost,randn(n),n)
    
H_ghost = diagm(-1 => β, 0 => α, 1 => β)
approx_eigvals_lanczos_ghost = sort(eigvals(H_ghost))

# Visualize true eigenvalues and Lanczos-generated eigenvalues
figure_ghost = Figure(resolution=(1000, 600))
ax = Axis(figure_ghost[1, 1], xlabel="Index", ylabel="Eigenvalue")
scatter!(ax, 1:n, true_eigvals_ghost, color=:blue, label = "True Eigenvalues")
scatter!(ax, 1:n, approx_eigvals_lanczos_ghost, color=:red, label = "Lanczos Eigenvalues")
axislegend(ax, position = :lt)
figure_ghost
save("ghost_eigenvalues.png", figure_ghost)
import Pkg.instantiate
instantiate()
using BenchmarkTools: @ballocated
using Random: randperm
using LinearAlgebra: I, norm
using CairoMakie
include("HW3_your_code.jl")

#----------------------------------------
# Problem a + b
#----------------------------------------
########################################
A = randn(20, 20) + 100 * I
b = randn(20)
reference_x = A \ b
unpivoted_LU!(A) 
substitution!(b, A)
@assert reference_x ≈ b

allocated_memory = @ballocated  unpivoted_LU!(A)
allocated_memory += @ballocated  substitution!(b, A)
@assert allocated_memory < 450

#----------------------------------------
# Problem c
# Generate families of random 𝑚 × 𝑚 matrices and vectors of length 𝑚. Plot as a function of the
# size 𝑚, the relative error of the solution obtained from your code in parts (a,b) and the growth
# factor introduced in problem 2. Report the floating point type used by your program.
#----------------------------------------
# YOUR CODE GOES HERE
size_list = 2:2:1000
error_list = zeros(size(size_list,1))
growth_factor_list = zeros(size(size_list,1))
for i = 1:length(size_list)
    m = size_list[i]
    A = randn(m, m) + 100 * I
    A_copy = copy(A)
    b = randn(m)
    b_copy = copy(b)
    A_max = maximum(abs.(A))
    reference_x = A \ b
    unpivoted_LU!(A) 
    substitution!(b, A)
    error_list[i] = norm(A_copy*b - b_copy) / norm(b_copy)
    U_max = 0
    for i = 1:m
        for j = 1:m
            if i < j
                if abs(A[i,j]) > U_max
                    U_max = abs(A[i,j])
                end
            end
        end
    end
    growth_factor_list[i] = U_max / A_max
end
fig, ax, pltobj = scatter(size_list, error_list, color = :red, label = "unpivoted LU")
ax.xlabel = "size of matrix (m)"
ax.ylabel = "relative error"
xlims!(minimum(size_list), maximum(size_list))
ax.title = "relative error of the solution obtained from unpivoted LU"
save("unpivoted_relative_error.pdf", fig)

fig, ax, pltobj = scatter(size_list, growth_factor_list, color = :red, label = "unpivoted LU")
ax.xlabel = "size of matrix (m)"
ax.ylabel = "growth factor"
xlims!(minimum(size_list), maximum(size_list))
ax.title = "growth factor of A from unpivoted LU"
save("unpivoted_growth_factor.pdf", fig)

#----------------------------------------
# Problem d
#----------------------------------------
########################################
A = (randn(4, 4) + 100 * I)[randperm(4), :]
b = randn(4)
reference_x = A \ b
P = pivoted_LU!(A) 
substitution!(b, A, P)
@assert reference_x ≈ b

#plotting
size_list = 2:2:1000
error_list = zeros(size(size_list,1))
growth_factor_list = zeros(size(size_list,1))
for i = 1:length(size_list)
    m = size_list[i]
    A = randn(m, m) + 100 * I
    A_copy = copy(A)
    b = randn(m)
    b_copy = copy(b)
    A_max = maximum(abs.(A))
    reference_x = A \ b
    P = pivoted_LU!(A) 
    substitution!(b, A, P)
    error_list[i] = norm(A_copy*b - b_copy) / norm(b_copy)
    U_max = 0
    for i = 1:m
        for j = 1:m
            if i < j
                if abs(A[i,j]) > U_max
                    U_max = abs(A[i,j])
                end
            end
        end
    end
    growth_factor_list[i] = U_max / A_max
end

fig, ax, pltobj = scatter(size_list, error_list, color = :red, label = "pivoted LU")
ax.xlabel = "size of matrix (m)"
ax.ylabel = "relative error"
xlims!(minimum(size_list), maximum(size_list))
ax.title = "relative error of the solution obtained from pivoted LU"
save("pivoted_relative_error.pdf", fig)

fig, ax, pltobj = scatter(size_list, growth_factor_list, color = :red, label = "pivoted LU")
ax.xlabel = "size of matrix (m)"
ax.ylabel = "growth factor"
xlims!(minimum(size_list), maximum(size_list))
ax.title = "growth factor of A from pivoted LU"
save("pivoted_growth_factor.pdf", fig)

#----------------------------------------
# Problem e
#----------------------------------------
# YOUR CODE GOES HERE
size_list = 2:2:1000
error_list = zeros(size(size_list,1))
error_list_julia = zeros(size(size_list,1))
for i = 1:length(size_list)
    m = size_list[i]
    A = growth_matrix(m)
    A_copy = copy(A)
    b = randn(m)
    b_copy = copy(b)
    reference_x = A \ b
    P = pivoted_LU!(A) 
    substitution!(b, A, P)
    error_list[i] = norm(A_copy*b - b_copy) / norm(b_copy)
    error_list_julia[i] = norm(A_copy*reference_x - b_copy) / norm(b_copy)
end

fig, ax, pltobj = scatter(size_list, error_list, color = :red, label = "Pivoted LU")
ax.xlabel = "size of matrix (m)" 
ax.ylabel = "relative error"
xlims!(minimum(size_list), maximum(size_list))
ax.title = "Pivoted LU vs Julia built-in operator(\\)"
scatter!(size_list, error_list_julia, color = :blue, label = "Built-in operator")
axislegend(ax)
save("highGrowth_both_relative_error.pdf", fig)




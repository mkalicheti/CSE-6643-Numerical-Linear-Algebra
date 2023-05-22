import Pkg.instantiate
instantiate()
using BenchmarkTools: @ballocated
using LinearAlgebra: I, norm, istriu, triu, qr
using CairoMakie
include("HW4_your_code.jl")


#----------------------------------------
# Problem a 
#----------------------------------------
########################################
A = randn(30, 20) 
b = randn(30)
Q, R = classical_gram_schmidt(A) 
@assert Q' * Q ≈ I
@assert Q * R ≈ A

#----------------------------------------
# Problem b 
#----------------------------------------
########################################
A = randn(30, 20) 
b = randn(30)
Q, R = modified_gram_schmidt(A) 
@assert Q' * Q ≈ I
@assert Q * R ≈ A

#----------------------------------------
# Problem c
#----------------------------------------
########################################
A = randn(25, 20) 
allocated_memory = @ballocated  householder_QR!(A)
@assert allocated_memory == 0
A = randn(25, 20)
A_copy = copy(A)
true_R = Matrix(qr(A).R)
householder_QR!(A)
# Checks if the R part of the factorization is correct
@assert vcat(true_R, zeros(5,20)) ≈ triu(A)

#----------------------------------------
# Problem d
#----------------------------------------
########################################
# Testing for memory allocation:

## Note: the div function is incomplete: I was not able to get it to work
A = randn(25, 20) 
householder_QR!(A)
QR = A
x = randn(20)
b = randn(25)
out_mul = randn(25)
out_div = randn(20)

allocated_memory_mul = @ballocated  householder_QR_mul!(out_mul, x, QR)
# allocated_memory_div = @ballocated  householder_QR_div!(out_div, b, QR)
@assert allocated_memory_mul == 0
# @assert allocated_memory_div == 0

# Testing for correctness:
A = randn(5, 2) 
A_copy = copy(A)
x = randn(2)
b = randn(5)
b_copy = copy(b)
out_mul = randn(5)
out_mul_copy = copy(out_mul)
out_div = randn(2)
out_div_copy = copy(out_div)
true_mul = A * x 
true_div = A \ b 

householder_QR!(A)
QR = A
householder_QR_mul!(out_mul, x, QR)
householder_QR_div!(out_div, b, QR)

# checks whether the results are approximately correct
@assert true_mul ≈ out_mul
#@assert true_div ≈ out_div


#----------------------------------------
# Problem e
#----------------------------------------
# YOUR CODE GOES HERE

size_list = 2:2:2^6
CGS_err_list = zeros(size(size_list,1))
MGS_err_list = zeros(size(size_list,1))
house_err_list = zeros(size(size_list,1))
#Computing Errors
for i = 1:size(size_list,1)
    m = size_list[i]
    #True R matrix
    A_orig = randn(m, m)
    true_R = qr(A_orig).R
    #Classical Gram-Schmidt
    A = copy(A_orig)
    Q, R = classical_gram_schmidt(A) 
    CGS_err_list[i] = maximum(abs.(true_R - R))
    #Modified Gram-Schmidt
    A = copy(A_orig)
    Q, R = modified_gram_schmidt(A) 
    MGS_err_list[i] = maximum(abs.(true_R - R))
    #Householder Reflection
    A = copy(A_orig)
    householder_QR!(A)
    house_err_list[i] = maximum(abs.(true_R - triu(A)))
end
#Plotting
fig, ax, pltobj = scatter(size_list, CGS_err_list, color = :green, label = "Classical Gram-Schmidt")
scatter!(size_list, MGS_err_list, color = :red, label = "Modified Gram-Schmidt")
scatter!(size_list, house_err_list, color = :blue, label = "Householder reflection")
ax.xlabel = "Size of matrix"
ax.ylabel = "Error in R matrix"
ax.title = "Classical Gram-Schmidt, Modified Gram-Schmidt, and Householder reflection : Comparison of Stability"
axislegend(ax)
save("Error_Plot.pdf", fig)

#Instability of Classical Gram-Schmidt -  explored through edge cases
m = 32

#Example 1: Hilbert Matrix
println("---Eg.1---")
hilb = m->[1//(BigInt(i)+j-1) for i in 1:m, j in 1:m]
A_og = 0.000001*I + hilb(m)

true_R = qr(A_og).R

A = copy(A_og)
Q, R = classical_gram_schmidt(A) 
CGS_err = norm(I - Q'*Q)
println("Classical Gram-Schmidt's Error in R: ", CGS_err)

A = copy(A_og)
Q, R = modified_gram_schmidt(A) 
MGS_err = norm(I - Q'*Q)
println("Modified Gram-Schmidt's Error in R: ", MGS_err)

A = copy(A_og)
householder_QR!(A)
house_err = maximum(abs.(true_R - triu(A)))
println("Householder Reflection's Error in R: ", house_err)

#Example 2: Linear Dependency
println("---Eg.2---")
A_og = zeros(m, m)
for i = 1:m
    A_og[1, i] = 1
end
for i = 2:m
    A_og[i, 1] = 1/1000
end
for i=2:m
    A_og[i, i] = 1/1000
end
A_og[1:m,m] = 3*A_og[1:m, 1] + 4*A_og[1:m, 2] # Make linearly dependent

true_R = qr(A_og).R

A = copy(A_og)
Q, R = classical_gram_schmidt(A) 
CGS_err = norm(I - Q'*Q)
println("Classical Gram-Schmidt's Error in R: ", CGS_err)

A = copy(A_og)
Q, R = modified_gram_schmidt(A) 
MGS_err = norm(I - Q'*Q)
println("Modified Gram-Schmidt's Error in R: ", MGS_err)

A = copy(A_og)
householder_QR!(A)
house_err = maximum(abs.(true_R - triu(A)))
println("Householder Reflection's Error in R: ", house_err)

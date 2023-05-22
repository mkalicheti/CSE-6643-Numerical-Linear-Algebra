import Pkg.instantiate
instantiate()
using BenchmarkTools: @btime, @belapsed, @ballocated
using LinearAlgebra: mul!
using CairoMakie
include("HW2_your_code.jl")


#----------------------------------------
# Problem a
#----------------------------------------
########################################
# Change size_list to choose matrix 
# sizes suitable for your system
size_list = 2 .^ (6 : 12)
########################################
time_kernel = Float64[]
time_BLAS = Float64[]
for sz in size_list
    global A = rand(sz, sz)
    global B = rand(sz, sz)
    global C = rand(sz, sz)
    if sz < 1000
        push!(time_kernel, @belapsed add_to_A_B_times_C!(A, B, C))
    else 
        push!(time_kernel, @elapsed add_to_A_B_times_C!(A, B, C))
    end
    push!(time_BLAS, @belapsed mul!(A, B, C))
end
pl_dense = scatter(size_list, time_kernel, label="custom kernel", axis=(yscale=log10, xscale=log10, xlabel="problem size", ylabel="time in seconds"))
scatter!(size_list, time_BLAS, label="BLAS")
axislegend(position=:lt)
save("performance_per_size.pdf", pl_dense)
#----------------------------------------


#----------------------------------------
# Problem b
#----------------------------------------
# Testing blocked matmul function
A = zeros(2014, 301)
B = randn(2014, 1037)
C = randn(1037, 301)
println("Testing for memory allocation of blocked matmul...")
allocated_memory = @ballocated add_to_A_B_times_C!(A, B, C, 301)
@assert allocated_memory == 0
println("No memory allocated, good!")
println("Testing correctness of blocked matmul...")
A .= 0; add_to_A_B_times_C!(A, B, C, 301)
@assert A ≈ B * C
println("Result correct, good!")
#----------------------------------------


#----------------------------------------
# Problem c
#----------------------------------------
########################################
# Change size_list to choose matrix 
# sizes suitable for your system
l = 6010
m = 6011
n = 6070
# Modify sizes suitable for your system
block_sizes = 2 .^ (2 : 12)
########################################
time_blocked = Float64[]
for bks in block_sizes
    A = zeros(l, n)
    B = randn(l, m)
    C = randn(m, n)
    push!(time_blocked, @elapsed add_to_A_B_times_C!(A, B, C, bks))
end
pl_blocked = scatter(block_sizes, time_blocked, label="blocked", axis=(xlabel="block size", ylabel="time", yscale=log10))
axislegend(position=:lt)
save("performance_blocked.pdf", pl_blocked)
#----------------------------------------

#----------------------------------------
# Problem d
#----------------------------------------
# Testing matmul function
A = zeros(2014, 301)
B = randn(2014, 1037)
C = randn(1037, 301)
println("Testing for memory allocation of oblivious matmul...")
allocated_memory = @ballocated oblivious_add_to_A_B_times_C!(A, B, C, 64)
@assert allocated_memory == 0
println("No memory allocated, good!")
println("Testing correctness of oblivious matmul...")
A .= 0; oblivious_add_to_A_B_times_C!(A, B, C, 64)
@assert A ≈ B * C
println("Result correct, good!")
#----------------------------------------

#----------------------------------------
# Problem e
#----------------------------------------
time_oblivious = Float64[]
for bks in block_sizes
    A = zeros(l, n)
    B = randn(l, m)
    C = randn(m, n)
    push!(time_oblivious, @elapsed oblivious_add_to_A_B_times_C!(A, B, C, bks))
end
scatter!(block_sizes, time_oblivious, label="oblivious")
axislegend(position=:lt)
save("performance_blocked_and_oblivious.pdf", pl_blocked)
#----------------------------------------
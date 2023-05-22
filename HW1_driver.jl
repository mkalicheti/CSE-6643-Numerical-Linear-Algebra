# DO NOT MODIFY THIS FILE!
using BenchmarkTools: @btime
include("HW1_your_code.jl")
# Testing matvec function
m = 1000
n = 1700
u = randn(m)
A = randn(m, n)
v = randn(n)

println("Timing the matvec function...")
@btime u_is_A_times_v!(u, A, v)
println("Testing the matvec function...")
@assert u ≈ A * v
println("Matvec function passed test!")

# Testing matmul function
l = 1100
m = 1000
n = 1700
A = randn(l, n)
B = randn(l, m)
C = randn(m, n)

println("Timing the matmul function...")
@btime A_is_B_times_C!(A, B, C)
println("Testing the matmul function...")
@assert A ≈ B * C
println("Matmul function passed test!")
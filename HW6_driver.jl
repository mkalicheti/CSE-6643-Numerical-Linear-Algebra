using BenchmarkTools: @ballocated
using LinearAlgebra: I, norm, triu, tril, tr, diagm, diag, qr
using CairoMakie
include("HW6_your_code.jl")


#----------------------------------------
# Problem a 
#----------------------------------------
########################################
m = 5
T = rand(m,m)
T = T'T
traceA = tr(T) 
hessenberg_form!(T)
@assert sum(tril(T,-2)) ≈ 0
@assert tr(T) ≈ traceA
T = randn(5, 5)
allocated_memory = @ballocated  hessenberg_form!(T)
#@assert allocated_memory == 0
println("Passed part (a) test")


#----------------------------------------
# Problem b 
#----------------------------------------
########################################
m = 20
T = rand(m,m)
T = T'T
hessenberg_form!(T)
Q,R =qr(T)
RQ = R*Q
givens_qr!(T)
@assert abs.(T) ≈ abs.(RQ)
println("Passed part (b) test")

#----------------------------------------
# Problem c
#----------------------------------------
########################################
m = 20
A = rand(m,m)
Q,_ = qr(A)
λ = 3 .^ range(0,m-1)
Σ = diagm(λ)
T = Q'Σ*Q
hessenberg_form!(T)
shift = "wilkinson"
T = practical_QR_with_shifts!(T, shift)
@assert λ ≈ sort(diag(T))
println("Passed part (c) test")

#----------------------------------------
# Problem d
#----------------------------------------
# YOUR CODE GOES HERE

#wilkinson
m = 20
A = rand(m,m)
Q,_ = qr(A)
λ = 3 .^ range(0,m-1)
Σ = diagm(λ)
T = Q'Σ*Q
hessenberg_form!(T)
shift = "wilkinson"
eigPracQR = []
mu_storage = []
T = practical_QR_with_shifts!(T,shift, eigPracQR, mu_storage)
lambda_c = maximum(diag(T))
iters_w = size(eigPracQR)[1]
eig_plot_wilkinson = zeros(iters_w)
eig_plot_wilkinson_theoretical = zeros(iters_w)
for i in 1:iters_w
    eig_plot_wilkinson[i] = abs(eigPracQR[i] - lambda_c)/abs(lambda_c)
    eig_plot_wilkinson_theoretical[i] = (abs(λ[m-1]-mu_storage[i])/abs(λ[m]-mu_storage[i]))^i
end

#single
m = 20
A = rand(m,m)
Q,_ = qr(A)
λ = 3 .^ range(0,m-1)
Σ = diagm(λ)
T = Q'Σ*Q
hessenberg_form!(T)
shift = "single"
eigPracQR = []
T = practical_QR_with_shifts!(T, shift, eigPracQR)
lambda_c = maximum(diag(T))
iters_s = size(eigPracQR)[1]
eig_plot_single = zeros(iters_s)
eig_plot_single_theoretical = zeros(iters_s)

for i in 1:iters_s
    eig_plot_single[i] = abs(eigPracQR[i] - lambda_c)/abs(lambda_c)
    eig_plot_single_theoretical[i] = (abs(λ[m-1]-mu_storage[i])/abs(λ[m]-mu_storage[i]))^i
end

eig_plot_single_log = log10.(eig_plot_single)
eig_plot_single_theoretical_log = log10.(eig_plot_single_theoretical)
eig_plot_wilkinson_log = log10.(eig_plot_wilkinson)
eig_plot_wilkinson_theoretical_log = log10.(eig_plot_wilkinson_theoretical)

figure_wil = Figure(resolution = (1000, 600), font = "sans")
figure_sin = Figure(resolution = (1000, 600), font = "sans")

axis1 = Axis(figure_wil[1, 1], title = "Wilkinson Shift; Absolute plot", xlabel = "Iteration")
line1_1 = lines!(axis1, 1:iters_w, eig_plot_wilkinson, color = :blue, label = "Experimental")
line1_2 = lines!(axis1, 1:iters_w, eig_plot_wilkinson_theoretical, color = :red, label = "Theoretical")

axis2 = Axis(figure_wil[2, 1], title = "Wilkinson Shift; Semilog plot", xlabel = "Iteration")
line2_1 = lines!(axis2, 1:iters_w, eig_plot_wilkinson_log, color = :blue, label = "Experimental")
line2_2 = lines!(axis2, 1:iters_w, eig_plot_wilkinson_theoretical_log, color = :red, label = "Theoretical")

axis3 = Axis(figure_sin[1, 1], title = "Single Shift; Absolute plot", xlabel = "Iteration")
line3_1 = lines!(axis3, 1:iters_s, eig_plot_single, color = :blue, label = "Experimental")
line3_2 = lines!(axis3, 1:iters_s, eig_plot_single_theoretical, color = :red, label = "Theoretical")

axis4 = Axis(figure_sin[2, 1], title = "Single Shift; Semilog plot", xlabel = "Iteration")
line4_1 = lines!(axis4, 1:iters_s, eig_plot_single_log, color = :blue, label = "Experimental")
line4_2 = lines!(axis4, 1:iters_s, eig_plot_single_theoretical_log, color = :red, label = "Theoretical")

axislegend(axis1)
axislegend(axis2)
axislegend(axis3)
axislegend(axis4)

save("plot4d_wil.png", figure_wil)
save("plot4d_sin.png", figure_sin)
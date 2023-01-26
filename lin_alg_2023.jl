using LinearAlgebra
using Plots
using Serialization
using Statistics
using LaTeXStrings

let # to establish local scope so variables don't be annoying

# Function definitions

	# poly_approx takes x sample points xs, y sample points ys, and an AbstractVector ns
	# that contains all of the numbers n for which x^n should be in the set of basis 
	# functions.
	#
	# returns a Tuple containing the expansion coefficients, the approximation as applied to
	# the sampled x values, and the RMS error of the approximation on the sampled points.
	function poly_approx(xs, ys::Vector, ns::AbstractVector)
		# generating basis functions
		polys = [x ^ j for x in xs, j in ns]

		# solving least-squares problem for coefficients
		coefs = polys \ ys

		# evaluate polynomial on the sample points
		approx_ys = [sum(x .^ Vector(ns) .* coefs) for x in xs]

		# calculate RMS error
		rmserr = sqrt(mean((approx_ys .- ys) .^ 2 ))

		return coefs, approx_ys, rmserr
	end

	# poly_approx can also take a function so that you don't have to evaluate it
	# beforehand on the sample points
	poly_approx(xs, f::Function, ns) = poly_approx(xs, f.(xs), ns)


# Problem 1: Approximating our first function, sin(x)

	xs = LinRange(0, pi/2, 1_000)

	_, sinapp_ys, sinapp_err = poly_approx(xs, sin, 0:7)

	@show sinapp_err # = 1.2159061878334831e-8


	#plot(xs, sin.(xs); label=L"sin(x)", color=:blue)
	#p1 = plot!(xs, sinapp_ys; label="Degree-5 approx.", color=:green)

	plot(xs, sin.(xs) .- sinapp_ys; label="Approx. error", color=:red)
	#plot(p1, p2, layout=@layout [a ; b])

	savefig("sinapp.pdf")

# Problem 2: Quantifying the relationship between error and polynomial degree

	max_app_errs = 0:20
	app_errs = [poly_approx(xs, sin, 0:max_app_err)[3] for max_app_err in max_app_errs]

	scatter(max_app_errs, app_errs; yaxis=:log, label="Approx. errors", xaxis="Polynomial deg.", title=L"$\sin(x)$ approximation error")
	plot!(max_app_errs, fill(1.e-10, 21); label="1e-10")
	plot!(max_app_errs, fill(1.e-20, 21); label="1e-20")

	savefig("sinapperrs.pdf")

# Problem 3: Antisymmetric polynomials

	anti_max_app_errs = 1:2:20
	antisym_app_errs = [poly_approx(xs, sin, 1:2:max_app_err)[3] for max_app_err in anti_max_app_errs]

	scatter(max_app_errs, app_errs; yaxis=:log, label="All degs.", xaxis="Polynomial deg.", title=L"$\sin(x)$ approximation error")
	scatter!(anti_max_app_errs, antisym_app_errs; label="Antisym.")
	plot!(max_app_errs, fill(1.e-10, 21); label="1e-10")
	plot!(max_app_errs, fill(1.e-20, 21); label="1e-20")

	savefig("sinapperrsanti.pdf")

# Problem 4: Derivatives as linear operators

	f_coefs, _, f_err = poly_approx(xs, sin, 0:13)

	# deriv_mat creates the derivative matrix up to polynomial order n
	deriv_mat(n::Integer)= [i == j-1 ? j : 0. for i=0.:n, j=0.:n]

	# f′_coefs is the derivative calculated manually from the coefficients of the approx. to f = sin(x)
	f′_coefs = deriv_mat(13) * f_coefs

	# g_coefs are the coefficients to the polynom. approx. of the analytic derivative of f, cos(x)
	g_coefs, _, g_err = poly_approx(xs, cos, 0:13)
	sin_cos_err = sqrt(mean((g_coefs .- f′_coefs) .^ 2)) 

	@show f_err       # = 2.309736279836864e-16
	@show g_err       # = 3.3214652486332893e-16
	@show sin_cos_err # = 1.6752669448773883e-8

end

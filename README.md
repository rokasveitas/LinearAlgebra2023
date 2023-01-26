# LinearAlgebra2023

This is a repository containing my work for Project 2 of the Perimeter Scholars International course on Numerical Methods. All code is contained in `lin_alg_2023.jl`.  Here I describe how each section of the code works.

## Function definitions

The core function here is `poly_approx`. It solves the least-squares problem via the `\` backslash operation between a set of basis polynomials and the values of the desired function.  It also calculates the value of the approximating polynomial on the sample points and the root-mean-square (RMS) error of the approximation

## Problem 1

Here we approximate the function $\sin(x)$ on the interval $[0, \pi/2]$ over 1000 evenly-spaced sample points with a degree-seven polynomial.  The RMS error is `1.21e-8`.  We can plot the approximation error along the function domain:

![sinapp.png](https://github.com/rokasveitas/LinearAlgebra2023/plots/sinapp.png)

The dominant component of this error is that of a ninth-order polynomial, as we would expect for the antisymmetric $\sin(x)$ when we stop at seven.

## Problem 2

Now we'd like to look at how much error we can get rid of by going to higher polynomial degree.  We can carry out this same procedure at every integer order under 20 and plot the resulting errors:

![sinapperrs.png](https://github.com/rokasveitas/LinearAlgebra2023/plots/sinapperrs.png)

The error floor is on the order of `eps(Float64)`, which is `2.220446049250313e-16`.

## Problem 3

Because $\sin(x)$ is an antisymmetric function, we really didn't need to be using all of the even-degree monomials.  We can plot the errors for solely antisymmetric approximations at the same time as the previous data:

![sinapperrsanti.png](https://github.com/rokasveitas/LinearAlgebra2023/plots/sinapperrsanti.png)

We can see that there is some real information being lost by only using antisymmetric polynomials, and perhaps we could quantify that with an information criterion, but the error is still quite low and we reach the same floor eventually.

## Problem 4

Now we compare the derivative that we get by knowing how derivatives act on our polynomial approximation, and the derivative that we get by polynomial-approximating the analytic derivative of $\sin(x)$, which is $\cos(x)$.

Our data is as follows:
- The RMS error of the approximation $P\[\sin\](x)$ is `2.3e-16`
- The RMS error of the approximation $P\[\cos\](x)$ is `3.3e-16`
- The RMS error between $P\[\cos\](x)$ and $\frac{d}{dx} P\[\sin\](x)$ is `1.6e-8`.

We certainly accrue error, so these operations don't exactly commute, but this is as good as we can reasonably expect with floating-point linear algebra.

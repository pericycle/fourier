# fourier
Python code for computing and visualizing Fourier series approximation.

This code was written to create animations to demonstrate Fourier series
approximation to a class of 2nd year Engineering students at the University of Zimbabwe, Harare.

The module contains the following files:

# plotAnimFourier(case, L, n, delay)
case is a string for the case number.
default: case = 'wiggly'
case options: 'even', 'odd', 'wiggly', 'sin_abs', 'sin', 'sin_decay', 'abs', 'cos13', 'parabola',
                'parajump', 'sin_ext', 'jump', and 'line'
L is half the period, e.g. L = pi for the interval (-pi, pi)
n is the number of Fourier terms, e.g. n = 5
delay is a time.sleep parameter for the animation in seconds. default: delay = 1/16 

# plotFourier(case, L, n)
plots the function and Fourier series approximation with n terms for the case example

# fourierCoeffs(case, L, n)
computes the Fourier coefficients of the function case up to order n. 
returns two vectors a,b containing the values of the Fourier coefficients
a = [a_0, a_1, ...., a_n] and b = [b_1, b_2,.....,b_n]
we use the convention that f(x) ~ a_0 + a_1 cos(x pi/L) + b_1 sin(x pi/L) + ...
where the coefficient a_0 = 1/(2L) integral f(x) dx 

# simpson(f, a, b, n)
computes the integral of f from a to b using Simpson's rule
use n (number od points) odd so that the number of intervals is even





# fourier
Python code for computing and visualizing Fourier series approximation. <br/>

This code was written to create animations to demonstrate Fourier series
approximation to a class of 2nd year Engineering students at the <br/>
University of Zimbabwe, Harare. <br/>

The module contains the following files: <br/>

# plotAnimFourier(case, L, n, delay) <br/>
case is a string for the case number. <br/>
default: case = 'wiggly' <br/>

case options: 'even', 'odd', 'wiggly', 'sin_abs', 'sin', 'sin_decay', 'abs', 'cos13', 'parabola', <br/>
                'parajump', 'sin_ext', 'jump', and 'line'  <br/>
L is half the period, e.g. L = pi for the interval (-pi, pi) <br/>
n is the number of Fourier terms, e.g. n = 5 <br/>
delay is a time.sleep parameter for the animation in seconds. default: delay = 1/16  <br/>

# plotFourier(case, L, n) <br/>
plots the function and Fourier series approximation with n terms for the case example <br/>

# fourierCoeffs(case, L, n) <br/>
computes the Fourier coefficients of the function case up to order n. <br/>
returns two vectors a,b containing the values of the Fourier coefficients <br/>
a = [a_0, a_1, ...., a_n] and b = [b_1, b_2,.....,b_n]  <br/>
we use the convention that f(x) ~ a_0 + a_1 cos(x pi/L) + b_1 sin(x pi/L) + ... <br/>
where the coefficient a_0 = 1/(2L) integral f(x) dx  <br/>

# simpson(f, a, b, n)   <br/>
computes the integral of f from a to b using Simpson's rule  <br/>
use n (number od points) odd so that the number of intervals is even  <br/>

# funcgen(case) <br/>
this function returns an example function f





function y = Horner_Newton(c, a, x, n)
%
% Evaluates a(1)+a(2)(x-x(1))+a(3)(x-x(1))(x-x(2))+\cdots
% a(n+1)(x-x(1))(x-x(2))\cdots(x-x(n))
%
% Input c the value at which the polynomial is to be evaluated
% n the degree of the polynomial
% a (n+1) by 1 array that contains the coefficients of the
% polynomial
% x the x-components of the data points.
%
% Output y the value of the polynomial at x=c
%
y = a(n+1);
for i=1:n
y = y*(c-x(n+1-i))+a(n+1-i);
end
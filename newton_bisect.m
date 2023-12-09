function [x,fx,nf,af,bf] = newton_bisect(fname,a,b,tolx,nmax)
%
% Finds the root r of f(x)=0 using the best features of Bisection
% and Newton's methods.
%
% Input
% fname  handle of the function f(x) whose zero r is being computed.
%        it is also used to compute f'(x) directly
% a,b    he end points of the interval where the root of f(x) = 0
% tolx   a positive termination criteria
% nmax   maximum number of function evaluations
%
%
% Output
% x      the approximate root of the equation f(x)=0
% fx     the value of F at x
% nf     number of times f(x) and f'(x) have been evaluated
% af,bf  the end points of the final bracketing interval
%
%

x_1 = a;
x_2 = b;

%initialize the number of iterations
nf = 0;
N = 0;

%define tolerance and max # of iterations if not stated in the arguments
if nargin == 3
    tolx = 10^(-10); nmax = 50;
end

while (abs((x_1)-(x_2)) > tolx) && (N<nmax)
    %calculate the next iterate using newton step
    X2 = x_1 - (fname(x_1)./str2double(diff(fname(x_1)))); nf = nf +1;
    if (X2 >= x_1) && (X2 <= x_2)
        %check where the root is bracketed
        if fname(x_1)*fname(X2) <= 0
            x_2 = X2;
        else 
            x_1 = X2;
        end
    else
        %take a BIsection step when X2 is not in the interval
        X2 = ((x_1)+(x_2))/2;
        %check where the root is bracketed
        if fname(x_1)*fname(X2) <= 0
            x_2 = X2;
        else 
            x_1 = X2;
        end
    %count number of iterations
    N = N + 1;
    end
end
x = X2;
fx = fname(x);
af = x_1; bf = x_2;
disp(fprintf("x = %18.16f nf = %4f af = %4.8f bf = %4.8f",x,nf,af,bf))
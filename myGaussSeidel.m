function [X,N] = myGaussSeidel(A,b,x1,nmax,tol)
% Gauss-Seidel Method
%
% Method to solve Ax = b, where A=[aij]1<=i,j<=n, A is an nxn invertible
% matrix
% 
% Input
% 
% A    = the nxn invertible matrix
% b    = the right hand side of the equation ie the solution of Ax
% x1   = the initial guess
% nmax = maximum number of iterations
% tol  = tolerance to stop the iteration whenever the relative error 
%        norm(xi - xi-1)/norm(xi) is smaller than tol
%
%
% Output
%
% x = is the solution that will give Ax = b
% n = number of iterations the method took to converge
%

% Initiliaze values
N = 0;
n = length(b);
X=x1;
tolcheck = inf;


while (tolcheck > tol) && (N<=nmax)
    x0 = X;
    for i=1:n
        %initialize the sum of Aij*x^k-1
        temp1 = 0;
        %from 1 to i-1
        for j=1:i-1
            %summation of Aij*x^k from 1 to i-1
            temp1 = temp1 + (A(i,j)*X(j));
        end
        temp = 0;
        %from i+1 to n
        for j=i+1:n
            % so we dont add Aii
            if j~=i
                %add the summation of Aij*x^k-1, j!=i
                temp = temp + (A(i,j)*x0(j));
            end
        end
        % b - summation of Aij*x^k-1 - summation of Aij*x^k, j!=i divided by Aii
        X(i) = (b(i) - temp - temp1)/A(i,i);
    end
    %check is x^k-1 is a zero vector, if so skip or else itll give us 1
    if any(x0 ~= 0) && (N ~= 0)
        tolcheck = (norm((X-x0),1)).*(1./norm(X,1))
    end
    %add number of iterations
    N = N+1;
end


function a = newton_dd(X,F,n)
% this is a function that will give the coeeficients of the NPI using the given data points, and
% give out a list of coefficients for Newton's polynomial interpolation (a)
% X is a list of x-values of the data points
% F is a list of y or f(x) values of the data points
% n is the degree of the polynomial
    n = n + 1;
    a(1)=F(1);
    for i = 1:n - 1
        trix(i,1) = ((F(i+1)-F(i)) / (X(i+1)-X(i)) );
        % first column of the newton poly interpolation "matrix"
        % where f[x0,x1] f[x1,x2]........f[xn-1,xn]
    end
    for j = 2:n - 1
        for i = 1:n - j
            trix(i,j) = (trix(i+1,j-1) - trix(i,j-1))/(X(j+i)-X(i));
        end
    end
    trix
    for j = 2:n
        a(j)= trix(1,j-1);
    end
end

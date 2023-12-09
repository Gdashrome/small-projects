function plot_newton_poly(c,X,F)
    n = numel(X) - 1;
    disp([X; F])
    xlabel('x')
    ylabel('y')

    a = newton_dd(X,F,n);
    for i = 1:c
        x2(i) = X(1) + (i * ((abs(X(1)) + abs(X(n+1)))*(1/c)));
    end
    %divides the range of X in 100 different xi
    for i = 1:numel(x2)
        y2(i) = Horner_Newton(x2(i), a, X, n);
    end
    %evaluates each x2i using the Horner Newton function with the newton
    %coefficients taken from newton_dd function
    
    plot(x2,y2,'r-')
    %plot the coordinates
end
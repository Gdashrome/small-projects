function sum = Composite_Simpson(fname,a,b,n)
% this is a function that applies Simpson Quadrature Rule over the
% subdivided interval [a,b]
%
% fname is the function whose definite integral is being computed
% a and b are the end points of the interval of integration
% n is the number of intervals into which [a,b] gets subdivided

h = (a-b)./n;

intervals_a(1,1) = a;
for i = 2:n
    intervals_a(i,1) = (a+(i-1).*h);
end
intervals_a(n+1,1) = b;
intervals_a
sum = 0;
for i = 1:(numel(intervals_a) - 1)
    v = ((((intervals_a(i+1)-intervals_a(i))./6).* fname(intervals_a(i))) + ((2.*(intervals_a(i+1)-intervals_a(i))./3).*fname((intervals_a(i)+intervals_a(i+1))/2)) + (((intervals_a(i+1)-intervals_a(i))/6).*fname(intervals_a(i+1))));
    sum = sum + v;
end
disp sum;

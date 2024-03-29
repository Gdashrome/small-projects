x = -5:0.01:5;
y = 1./(1.+ x.^2);

X = [-5;0;5];
F = [1/26;1;1/26];
c = 100;

X1 = [-5;-2.5;0;2.5;5];
F1 = [1/26;4/29;1;4/29;1/26];

X2 = [-5;-10/3;-5/3;0;5/3;10/3;5];
F2 = [1/26;9/109;9/34;1;9/34;9/109;1/26];

X3 = [-5;-15/4;-5/2;-5/4;0;5/4;5/2;15/4;5];
F3 = [1/26;16/241;4/29;16/41;1;16/41;4/29;16/241;1/26];

hold on
plot(x,y,"b")
plot_newton_poly(c,X,F)
plot_newton_poly(c,X1,F1)
plot_newton_poly(c,X2,F2)
plot_newton_poly(c,X3,F3)
hold off
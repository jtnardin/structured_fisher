%RD_sim written 10-10-16 by JTN to simulate RD equation given parameters.

function w = RD_sim(D,lambda,t,x,IC)


    f = @(u) lambda*u.*(1-u);

    xn = length(x);
    dx = x(2)-x(1);
    x_int = 2:xn-1;

    A = D/dx^2*sparse([x_int x_int x_int],[x_int-1 x_int x_int+1],[ones(1,xn-2) -2*ones(1,xn-2) ones(1,xn-2)],xn,xn);
    A_bd = D/dx^2*sparse([1 1 xn xn],[1 2 xn-1 xn],[-2 2 2 -2],xn,xn);

    A = A+A_bd;

    [t,w] = ode45(@(t,w) RD_pde(t,w,f,A),t,IC);

end
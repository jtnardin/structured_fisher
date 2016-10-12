%RD_sim written 10-12-16 by JTN to simulate RD equation given parameters.

function w = RD_sim_nonaut(D1,D2,l1,l2,t,x,IC)

    %ensure psi varies between 0 and 0.3 (initial support)
    psi = @(t) min(max(0,1./(1+exp(2-2*cos(t)))),0.3);

    D = @(t) (D2 + (D1 - D2)*(10/3)*psi(t));

    lambda = @(t) l2 + (l1 - l2)*(10/3)*psi(t);
    
    f = @(u,t) lambda(t)*u.*(1-u);

    xn = length(x);
    dx = x(2)-x(1);
    x_int = 2:xn-1;

    A = @(t) D(t)/dx^2*sparse([x_int x_int x_int],[x_int-1 x_int x_int+1],[ones(1,xn-2) -2*ones(1,xn-2) ones(1,xn-2)],xn,xn);
    A_bd = @(t) D(t)/dx^2*sparse([1 1 xn xn],[1 2 xn-1 xn],[-2 2 2 -2],xn,xn);

    A = @(t) A(t) + A_bd(t);

    [t,w] = ode45(@(t,w) RD_pde_nonaut(t,w,f,A),t,IC);

end
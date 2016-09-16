%fisher_activity_structure_no_diff.m written 5-30-16 by JTN to simulate
%u_t + (g(m)u)_m = u(1-u). For now, g(m) = exp(-m), but subject to change
%in future.

n = 201;
dt = 1e-3;

x = linspace(0,5,n);
dx = x(2) - x(1);
t = 0:dt:10;

xn = length(x);
tn = length(t);

x_int = 2:xn-1;

D = 0;

% V = exp(-x(x_int)); %g(m) = exp(-m)
V = x(x_int); %g(m) = m
% V = x(x_int).^(1/2); %g(m) = m


theta = 0.5;

Dc = D*dt/dx^2;
Vc = V'*dt/dx;

u0 = 1;

%initial condition
% IC = @(x) u0*exp(-x);
% IC = @(x) 0.9*(x==0);
IC = @(x) u0*(x<=1);

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation

A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc+-v+v.*sw/2); ...
    (2*Dc+v-v.*se/2-v.*sw/2); (-Dc+v.*se/2)],xn,xn);

%velocity always positive here
%A_neg = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc-v.*sw/2); ...
    %(2*Dc+v.*se/2+v.*sw/2-v); (-Dc+v-v.*se/2)],xn,xn);

    
%initialize
u = zeros(xn,tn);
u(:,1) = IC(x);

tic

for i = 2:tn
    
    r_e = (u(x_int,i-1) - u(x_int-1,i-1))./(u(x_int+1,i-1) - u(x_int,i-1));    
    r_w = (u(x_int(2:end)-1,i-1) - u(x_int(2:end)-2,i-1))./(u(x_int(2:end),i-1) - u(x_int(2:end)-1,i-1));
    r_w = [-1;r_w];
    
    %eliminate NaN values (0/0 -- not steep!)
    r_e(isnan(r_e)) = 1;
    r_w(isnan(r_w)) = 1;
    
    %set inf values to large value
    r_e(isinf(r_e)) = 100;
    r_w(isinf(r_w)) = 100;

    %compute interior points
    u(:,i) = (speye(xn) + theta*A_pos(sigma(r_e),sigma(r_w),Vc,x_int,1))\((speye(xn) - (1-theta)*A_pos(sigma(r_e),sigma(r_w),Vc,x_int,1))*u(:,i-1) + dt*u(:,i-1).*(1-u(:,i-1)));
   
    
%     %boundary conditions
     %u(1,i) = u(1,i-1) + dt*(2*u(1,i-1)-u(1,i-1).^2); %u_m(m=0)=0
     u(1,i) = u(1,i-1) + dt/dx*(u(1,i-1)-u(2,i-1))+dt*(u(1,i-1)*(1-u(1,i-1))); %u_t + u_m = u(1-u)

end

toc

figure
hold on

colors = 'bgrkm';

j=1;
for i = 1:250:1250
    plot(x,u(:,i),colors(j))
    j=j+1;
end

xlabel('m')
ylabel('u(t,m)')
title('Numerical simulations of (1) over time')

j=1;
for i = 1:250:1250
    leg_times(j)=t(i);
    j=j+1;
end

legend(num2str(leg_times'),'location','southwest')

% exportfig(gcf,['u_num_u0_09_neumann.eps'])

figure

contour(t,x,u,11)
colorbar
hold on
plot(t,log(1+t),'--')
plot(t,log(exp(1)+t),'k--')

xlabel('t')
ylabel('m')
title('Contours of numerical simulations for (1)')

h=legend('contours','$\sigma ^{-1}(t,0)$','$\sigma ^{-1}(t,1)$','location','northeast');
set(h,'interpreter','latex')

% exportfig(gcf,['u_num_contour_u0_09_neumann.eps'])

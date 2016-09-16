%fisher_diffusion_act.m written 9-6-16 by JTN to simulate
%u_t + (g(m)u)_m = u_xx u. For now, g(m) = exp(-m), but subject to change
%in future.

%need to better approximate the velocity -- how to account for places where
% g=0?

clear all; clc

%Construct vectors of independent variables
mn = 21; %number of m points
xn = 101; %number of x points
total = mn*xn;

dt = 1e-3; %time step
t = 0:dt:5;

m = linspace(0,1.1,mn);
dm = m(2) - m(1);

x = linspace(0,10,xn);
dx = x(2) - x(1);

[X,M] = meshgrid(x,m);

tn = length(t);

%construct boundary points, interior
xm_int = 1:mn*xn;

x_bd_0 = 1:mn;
x_bd_l = (xn-1)*mn+1:xn*mn;

x_bd = union(x_bd_0,x_bd_l);

m_bd_0 = 1:mn:mn*xn;
m_bd_1 = mn:mn:xn*mn;

m_bd = union(m_bd_0,m_bd_1);

bd = union(m_bd,x_bd);

xm_int(bd) = [];

%rate of cellular diffusion
D = 1;
%rate of MAPK activation
%V = .05*ones(length(m)-2,1);
V = (1 - m(2:end-1))/20; %g(m) = 1 - m

%measure of implicitness (1 = backward euler, 0 = forward euler, 1/2 =
%crank nicholson)
theta = 0.5;

%computational diffusion, velocity
Dc = D*dt/dx^2;
Vc = V'*dt/dm;
Vc = repmat(Vc,1,xn-2);
Vc = Vc(:);

%initial condition
u0 = 1;
IC = u0*(X<=6).*(X>=2).*(M>=.1).*(M<=.2);
% IC = u0*(X<=6).*(X>=2).*M.*exp(-M);
IC = IC(:);

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation

% A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc+-v+v.*sw/2); ...
%     (2*Dc+v-v.*se/2-v.*sw/2); (-Dc+v.*se/2)],total,total);

A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-v+v.*sw/2); ...
    (+v-v.*se/2-v.*sw/2); (v.*se/2)],total,total);

Dmat = Dc*sparse([xm_int xm_int xm_int],[xm_int-mn xm_int+mn xm_int],[ones(1,2*length(xm_int))...
    -2*ones(1,length(xm_int))],total,total);

Dmat_bd = + Dc*sparse([x_bd_0 x_bd_0 x_bd_l x_bd_l],[x_bd_0 x_bd_0+mn x_bd_l x_bd_l-mn],...
        2*[-ones(1,length(x_bd_0)) ones(1,length(x_bd_0)) -ones(1,length(x_bd_0))...
        ones(1,length(x_bd_0))],total,total);

Dmat = Dmat + Dmat_bd;
    
%velocity always positive here
%A_neg = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc-v.*sw/2); ...
    %(2*Dc+v.*se/2+v.*sw/2-v); (-Dc+v-v.*se/2)],xn,xn);

    
%initialize
u = zeros(total,tn);
u(:,1) = IC;

tic

for i = 2:tn
    
    r_e = (u(xm_int,i-1) - u(xm_int-1,i-1))./(u(xm_int+1,i-1) - u(xm_int,i-1));    
    r_w = (u(xm_int(2:end)-1,i-1) - u(xm_int(2:end)-2,i-1))./(u(xm_int(2:end),i-1) - u(xm_int(2:end)-1,i-1));
    r_w = [-1;r_w];
    
    %eliminate NaN values (0/0 -- not steep!)
    r_e(isnan(r_e)) = 1;
    r_w(isnan(r_w)) = 1;
    
    %set inf values to large value
    r_e(isinf(r_e)) = 100;
    r_w(isinf(r_w)) = 100;

    %compute interior points
    u(:,i) = (speye(total) + theta*A_pos(sigma(r_e),sigma(r_w),Vc,xm_int,1))\((speye(total)...
        - (1-theta)*A_pos(sigma(r_e),sigma(r_w),Vc,xm_int,1) + Dmat)*u(:,i-1)); %+ ...
        %dt*u(:,i-1).*(1-u(:,i-1)));
   
    
%     %boundary conditions
     %u(1,i) = u(1,i-1) + dt*(2*u(1,i-1)-u(1,i-1).^2); %u_m(m=0)=0
     %u(1,i) = u(1,i-1) + dt/dx*(u(1,i-1)-u(2,i-1))+dt*(u(1,i-1)*(1-u(1,i-1))); %u_t + u_m = u(1-u)

end

toc

%y for visualization
y = zeros(mn,xn,tn);
z = zeros(xn,tn);

for i = 1:tn
    y(:,:,i) = reshape(u(:,i),mn,xn);
    z(:,i) = dm*trapz(y(:,:,i));
end
figure
for i = 1:100:tn
    subplot(1,2,1)
    surf(x,m,y(:,:,i),'edgecolor','none')    
    title(['t = ' num2str(t(i))])
    axis([0 10 0 1 0 1])
    caxis([0,1])
    colorbar
    view(2)
    
    subplot(1,2,2)
    plot(x,z(:,i))    
    axis([0 10 0 1])
    pause(.125)
end

figure

plot(sum(y(:,:,ceil(end)),2))
hold on
plot(sum(y(:,:,ceil(end/2)),2))
plot(sum(y(:,:,1),2))

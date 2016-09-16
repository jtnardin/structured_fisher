%fisher_diffusion_act_incorp_ve_ve.m written 9-7-16 by JTN to simulate
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

m = linspace(0,1,mn);
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

%define activation modulus
g = @(m) (1-m)/10;
%must be a better way... but should work for now
for i = 1:mn-1
    V(i) = g((m(i)+m(i+1))/2);
end

%construct vectors for east and west velocity.
Ve = V(2:end);
Vw = V(1:end-1);
Vw_m1 = V(end);

%measure of implicitness (1 = backward euler, 0 = forward euler, 1/2 =
%crank nicholson)
theta = 0.5;

%computational diffusion, east velocity and west velocity
Dc = D*dt/dx^2;

Vec = Ve'*dt/dm;
Vec = repmat(Vec,1,xn-2);
Vec = Vec(:);

Vwc = Vw'*dt/dm;
Vwc = repmat(Vwc,1,xn-2);
Vwc = Vwc(:);

Vwc_m1 = Vw_m1'*dt/dm;
Vwc_m1 = repmat(Vwc_m1,1,xn-2);
Vwc_m1 = Vwc_m1(:);




%initial condition
u0 = 1;
% IC = u0*(X<=6).*(X>=2).*(M>=.2).*(M<=.4);
IC = u0*(X<=6).*(M<=0.5).*(M>=0.05).*exp(-M);
IC = IC(:);

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation

% A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc+-v+v.*sw/2); ...
%     (2*Dc+v-v.*se/2-v.*sw/2); (-Dc+v.*se/2)],total,total);

A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[vw.*(-1+sw/2); ...
    (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],total,total);


A_pos_m1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
    (-vw.*sw/2)],total,total);


Dmat = Dc*sparse([xm_int xm_int xm_int],[xm_int-mn xm_int+mn xm_int],[ones(1,2*length(xm_int))...
    -2*ones(1,length(xm_int))],total,total);

Dmat_bd = + Dc*sparse([x_bd_0 x_bd_0 x_bd_l x_bd_l],[x_bd_0 x_bd_0+mn x_bd_l x_bd_l-mn],...
        2*[-ones(1,length(x_bd_0)) ones(1,length(x_bd_0)) -ones(1,length(x_bd_0))...
        ones(1,length(x_bd_0))],total,total);

Dmat = Dmat + Dmat_bd;
    
A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-vw.*sw/2); ...
    (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],xn,xn);

    
%initialize
u = zeros(total,tn);
u(:,1) = IC;

tic

for i = 2:tn
    
%     r_e = (u(xm_int,i-1) - u(xm_int-1,i-1))./(u(xm_int+1,i-1) - u(xm_int,i-1));    
%     r_w = (u(xm_int(2:end)-1,i-1) - u(xm_int(2:end)-2,i-1))./(u(xm_int(2:end),i-1) - u(xm_int(2:end)-1,i-1));
%     r_w = [-1;r_w];
%     %boundary m = 1
%     r_w_m1 = (u(m_bd_1(2:end-1)-1,i-1) - u(m_bd_1(2:end-1)-2,i-1))./(u(m_bd_1(2:end-1),i-1) - u(m_bd_1(2:end-1),i-1));
%     
%     %eliminate NaN values (0/0 -- not steep!)
%     r_e(isnan(r_e)) = 1;
%     r_w(isnan(r_w)) = 1;
%     r_w_m1(isnan(r_w_m1)) = 1;
%     
%     %set inf values to large value
%     r_e(isinf(r_e)) = 100;
%     r_w(isinf(r_w)) = 100;
%     r_w_m1(isinf(r_w_m1)) = 100;

    [r_e,r_w,r_w_m1] = positive_sensor(u(:,i-1),xm_int,m_bd_1);
    
    
    %compute interior points
    u(:,i) = (speye(total) + theta*(A_pos(sigma(r_e),sigma(r_w),Vec,Vwc,xm_int,1)...
        + A_pos_m1(sigma(r_w_m1),Vwc_m1,m_bd_1(2:end-1),1)))\((speye(total)...
        - (1-theta)*(A_pos(sigma(r_e),sigma(r_w),Vec,Vwc,xm_int,1)...
        + A_pos_m1(sigma(r_w_m1),Vwc_m1,m_bd_1(2:end-1),1)) + Dmat)*u(:,i-1)); %+ ...
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
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:100:tn
    subplot(1,2,1)
    surf(x,m,y(:,:,i),'edgecolor','none')    
    title(['t = ' num2str(t(i))])
    axis([0 10 0 1 0 1])
    caxis([0,1])
    colorbar
    view([-1 1 1])
    
    subplot(1,2,2)
    plot(x,z(:,i))    
    axis([0 10 0 1])
    pause(.125)
end
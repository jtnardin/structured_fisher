%fisher_diffusion_act_incorp_ve_ve.m written 9-7-16 by JTN to simulate
%u_t + (g(m)u)_m = u_xx u. For now, g(m) = exp(-m), but subject to change
%in future.

%need to better approximate the velocity -- how to account for places where
% g=0?

% clear all; clc

%Construct vectors of independent variables
mn = 21; %number of m points
xn = 201; %number of x points
total = mn*xn;
dt = 1e-3; %time step
t = 0:dt:15;
m = linspace(0,1,mn);
dm = m(2) - m(1);
x = linspace(0,20,xn);
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
D = 0;

%rate of proliferation
lambda = 0;

%define activation modulus, signal factor, and response to signal factor
s = @(t) 1+sin(t);

f = @(s) s-1;

g = @(m) m.*(1-m);%(1-m)/4;

%now define the velocity (in the m-direction) given the above variables
V = zeros(mn-1,1);
for i = 1:mn-1
    V(i) = g((m(i)+m(i+1))/2);
end

%construct vectors for east and west velocity.
Ve = V(2:end);
Vw = V(1:end-1);
Vw_m1 = V(end);
Ve_m0 = V(1);

%measure of implicitness (1 = backward euler, 0 = forward euler, 1/2 =
%crank-nicholson)
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

Vec_m0 = Ve_m0'*dt/dm;
Vec_m0 = repmat(Vec_m0,1,xn-2);
Vec_m0 = Vec_m0(:);


% 
Vec = @(t) f(s(t))*Vec;
Vwc = @(t) f(s(t))*Vwc;
Vwc_m1 = @(t) f(s(t))*Vwc_m1;
Vec_m0 = @(t) f(s(t))*Vec_m0;



%initial condition
u0 = 1;
% IC = u0*(X<=6).*(X>=2).*(M>=.2).*(M<=.4);
IC = u0*(X<=10).*(X>=5).*(M>=0.8).*exp(-M);
IC = IC(:);



%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%Define sparse matrix as a function for computation
%need to define for when velocity is both positive and negative.

A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[vw.*(-1+sw/2); ...
    (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],total,total);

A_pos_m1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
    (-vw.*sw/2)],total,total);

A_pos_m0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*(1-1*se/2); ve.*se/2],total,total);

Dmat = Dc*sparse([xm_int xm_int xm_int],[xm_int-mn xm_int+mn xm_int],[ones(1,2*length(xm_int))...
    -2*ones(1,length(xm_int))],total,total);

Dmat_bd = + Dc*sparse([x_bd_0 x_bd_0 x_bd_l x_bd_l],[x_bd_0 x_bd_0+mn x_bd_l x_bd_l-mn],...
        2*[-ones(1,length(x_bd_0)) ones(1,length(x_bd_0)) -ones(1,length(x_bd_0))...
        ones(1,length(x_bd_0))],total,total);

Dmat = Dmat + Dmat_bd; %incorporate both interior and boundary for simplicity
    
A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-vw.*sw/2); ...
    (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],total,total);

A_neg_m1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[(-vw.*sw/2); ...
    (vw.*sw/2-vw)],total,total);

A_neg_m0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*se/2; (ve-ve.*se/2)],total,total);


%vector for integrating over m
integ_mat = [];
for i = 1:xn
    integ_mat = [integ_mat; i*ones(mn,1)];
end
%initialize
u = zeros(total,tn);
u(:,1) = IC;

tic

for i = 2:tn

    
    %Change sensors, A matrix based on velocity (whose sign depends on
    %f(s(t))
    if f(s(t(i))) >= 0 
        %sensors
        [r_e,r_w,r_w_m1,r_e_m1] = positive_sensor(u(:,i-1),xm_int,m_bd_1,m_bd_0);
        %construct matrices
        A_int = A_pos(sigma(r_e),sigma(r_w),Vec(t(i)),Vwc(t(i)),xm_int,1);
        A_m1 = A_pos_m1(sigma(r_w_m1),Vwc_m1(t(i)),m_bd_1(2:end-1),1);
        A_m0 = A_pos_m0(sigma(r_e_m1),Vec_m0(t(i)),m_bd_0(2:end-1),1);
    else
        %sensors
        [r_e,r_w,r_w_m1,r_e_m1] = negative_sensor(u(:,i-1),xm_int,m_bd_1,m_bd_0);
        %construct matrices
        A_int = A_neg(sigma(r_e),sigma(r_w),Vec(t(i)),Vwc(t(i)),xm_int,1);
        A_m1 = A_neg_m1(sigma(r_w_m1),Vwc_m1(t(i)),m_bd_1(2:end-1),1);
        A_m0 = A_pos_m0(sigma(r_e_m1),Vec_m0(t(i)),m_bd_0(2:end-1),1);
    end
    
    %integrate over m (just riemann for now)
    w = dm*accumarray(integ_mat,u(:,i-1));
    w = repmat(w',mn,1);
    w = w(:);
    
    %main computation
    u(:,i) = (speye(total) + theta*(A_int + A_m1 + A_m0))\((speye(total)...
        - (1-theta)*(A_int + A_m1 + A_m0) + Dmat)*u(:,i-1) + dt*lambda*u(:,i-1).*(0.2-w)); 
  
end

toc

%y for visualization
y = zeros(mn,xn,tn);
z = zeros(xn,tn);

for i = 1:tn
    y(:,:,i) = reshape(u(:,i),mn,xn);
    z(:,i) = dm*sum(y(:,:,i));
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:100:tn
    subplot(1,2,1)
    surf(x,m,y(:,:,i),'edgecolor','none')    
    title(['t = ' num2str(t(i)) ', f(s) = ' num2str(f(s(t(i))))])
    axis([0 x(end) 0 1 0 1])
    caxis([0,1])
    colorbar
    view([-1 1 1])
    
    subplot(1,2,2)
    plot(x,z(:,i))    
    axis([0 x(end) 0 1])
    pause(.125)
end
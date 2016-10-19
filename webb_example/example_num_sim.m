

%need to better approximate the velocity -- how to account for places where
% g=0?

clear all; clc

%Construct vectors of independent variables
mn = 40; %number of m points
total = mn;

tn = 300;

t = linspace(0,1,tn);
dt = t(2) - t(1);
tn = length(t);


m = linspace(1,3,mn);
dm = m(2) - m(1);

[T,M] = meshgrid(t,m);

%define boundary points, interior
m_int = 2:mn-1;

m_bd_0 = 1;
m_bd_1 = mn;

m_bd = union(m_bd_0,m_bd_1);

%define activation modulus
g = @(m) sqrt(m-1).*(3-m);

%get east and west velocity points -- probably a better way, but doing this
%for now.
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

Vec = Ve'*dt/dm;

Vwc = Vw'*dt/dm;

Vwc_m1 = Vw_m1'*dt/dm;

%initial condition
% IC = u0*(X<=6).*(X>=2).*(M>=.2).*(M<=.4);
IC = 10*(m-1.1).^2.*(1.4-m).^2.*(m>=1.1).*(m<=1.4) + (m-2).^2.*(2.5-m).^2.*(m>=2).*(m<=2.5);

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation

% A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc+-v+v.*sw/2); ...
%     (2*Dc+v-v.*se/2-v.*sw/2); (-Dc+v.*se/2)],total,total);

A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[vw.*(-1+sw/2); ...
    (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],total,total);


A_pos_m1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
    (-vw.*sw/2)],total,total);
    
%initialize
u = zeros(total,tn);
u(:,1) = IC;

tic

for i = 2:tn
    
    r_e = (u(m_int,i-1) - u(m_int-1,i-1))./(u(m_int+1,i-1) - u(m_int,i-1));    
    r_w = (u(m_int(2:end)-1,i-1) - u(m_int(2:end)-2,i-1))./(u(m_int(2:end),i-1) - u(m_int(2:end)-1,i-1));
    r_w = [-1;r_w];
    %boundary m = 1
    r_w_m1 = (u(m_bd_1-1,i-1) - u(m_bd_1-2,i-1))./(u(m_bd_1,i-1) - u(m_bd_1-1,i-1));
    
    %eliminate NaN values (0/0 -- not steep!)
    r_e(isnan(r_e)) = 1;
    r_w(isnan(r_w)) = 1;
    r_w_m1(isnan(r_w_m1)) = 1;
    
    %set inf values to large value
    r_e(isinf(r_e)) = 100;
    r_w(isinf(r_w)) = 100;
    r_w_m1(isinf(r_w_m1)) = 100;
    
    
    %compute interior points
    u(:,i) = (speye(total) + theta*(A_pos(sigma(r_e),sigma(r_w),Vec,Vwc,m_int,1)...
        + A_pos_m1(sigma(r_w_m1),Vwc_m1,m_bd_1,1)))\((speye(total)...
        - (1-theta)*(A_pos(sigma(r_e),sigma(r_w),Vec,Vwc,m_int,1)...
        + A_pos_m1(sigma(r_w_m1),Vwc_m1,m_bd_1,1)))*u(:,i-1)); %+ ...
        %dt*u(:,i-1).*(1-u(:,i-1)));
   
    
%     %boundary conditions
     %u(1,i) = u(1,i-1) + dt*(2*u(1,i-1)-u(1,i-1).^2); %u_m(m=0)=0
     %u(1,i) = u(1,i-1) + dt/dx*(u(1,i-1)-u(2,i-1))+dt*(u(1,i-1)*(1-u(1,i-1))); %u_t + u_m = u(1-u)

end

toc

%y for visualization
y = zeros(mn,tn);
y = reshape(u,mn,tn);

figure

surf(t,m,y,'edgecolor','none')

xlabel('t')
ylabel('s')

caxis([0,0.01])
colorbar

axis([0 1 1 3 0 0.013])

view(2)
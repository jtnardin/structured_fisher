%%%%% compare_analytic_numerical.m written 10-19-16 by JTN to compare
%%%%% analytical and numerical simulations of u_t + (g(m)u)_m = 0 for u =
%%%%% u(t,m)

clear all; clc

%%%%% First define things pertinent to both analytic and numerical
%%%%% solutions

%grid
t = linspace(0,15,3000);
s = linspace(0,1,40);
[T,S] = meshgrid(t,s);

%load in g, sigma, sigma_inv, etc. functions
[g,sigma,sigma_inv,s_t,f,int_f_s] = g_sigma_h_example1;


%initial condition
phi = IC_uniform(.05,.35);



%%%%% generate analytic solution based on (Webb 2008)
Soln = @(t,s) g(sigma_inv(-t,s))./(g(s)).*phi(sigma_inv(-t,s));%.*(sigma_inv(t,1)<s);
%Exact soln matrix
A = Soln(T,S);



%%%% Now compute the solution numerically via flux limiters

%define boundary points, interior
sn = length(s);
tn = length(t);
ds = s(2) - s(1);
dt = t(2)-t(1);

s_int = 2:sn-1;

s_bd_0 = 1;
s_bd_n = sn;
s_bd_1_int = 1;
s_bd_nm1_int = sn-2;


s_bd = union(s_bd_0,s_bd_n);


%get east and west velocity points
for i = 1:sn-1
    V(i) = g((s(i)+s(i+1))/2);
end

%construct vectors for east and west velocity.
Ve = V(2:end);
Vw = V(1:end-1);
Ve_s0 = V(1);
Vw_s1 = V(end);

Vec = @(t) Ve'*dt/ds;

Vwc = @(t) Vw'*dt/ds;

Vec_s0 = @(t) Ve_s0'*dt/ds;

Vwc_s1 = @(t) Vw_s1'*dt/ds;


%measure of implicitness (1 = backward euler, 0 = forward euler, 1/2 =
%crank nicholson)
theta = 0.5;


%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));


A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[vw.*(-1+sw/2); ...
    (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],sn,sn);

A_pos_s0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*(1-1*se/2); ve.*se/2],sn,sn);

A_pos_s1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
    (-vw.*sw/2)],sn,sn);
    
%initialize
u = zeros(sn,tn);
u(:,1) = phi(s);

tic

for i = 2:tn
    
    [r_e,r_w,r_w_s1,r_e_s1] = positive_sensor(u(:,i-1),s_int,s_bd_n,s_bd_0,s_bd_1_int,s_bd_nm1_int);
        %construct matrices
        A_int = A_pos(sigma(r_e),sigma(r_w),Vec(t(i)),Vwc(t(i)),s_int,1);
        A_s1 = A_pos_s1(sigma(r_w_s1),Vwc_s1(t(i)),s_bd_n,1);
        A_s0 = A_pos_s0(sigma(r_e_s1),Vec_s0(t(i)),s_bd_0,1);
    
%     
%     %compute interior points
%     u(:,i) = (speye(total) + theta*(A_pos(sigma(r_e),sigma(r_w),Vec,Vwc,s_int,1)...
%         + A_pos_m1(sigma(r_w_m1),Vwc_m1,m_bd_1,1)))\((speye(total)...
%         - (1-theta)*(A_pos(sigma(r_e),sigma(r_w),Vec,Vwc,s_int,1)...
%         + A_pos_m1(sigma(r_w_m1),Vwc_m1,m_bd_1,1)))*u(:,i-1)); %+ ...
%         %dt*u(:,i-1).*(1-u(:,i-1)));
   
    u(:,i) = (speye(sn) + theta*(A_int + A_s1 + A_s0))\(speye(sn)...
        - (1-theta)*(A_int + A_s1 + A_s0))*u(:,i-1); 
  
        

end

toc

%y numerical simulation aoln Matrix
y = reshape(u,sn,tn);



%%%%%%% plot errors between numeric and analytic


figure

surf(t,s,A-y,'edgecolor','none')

xlabel('t')
ylabel('s')

% caxis([0,0.01])
% colorbar

% axis([0 1 1 3 0 0.013])

view(2)

step = floor(tn/4);
count = 1;
figure
for i = step:step:tn
    subplot(2,2,count)
    
    plot(s,A(:,i),'b')
    hold on
    plot(s,y(:,i),'r')
    title(['t = ' num2str(t(i))])
    if i == step
        legend('analyt','num','location','northwest')
    end
    
    count = count + 1;
end

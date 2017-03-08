%fisher_paper_nonaut_ex3.m written 11-14-16 by JTN to simulate
%nonautonomous Fisher's Equation for example 3


clear all; clc


%Construct vectors of independent variables
mn = 81; %number of m points
xn = 301; %number of x points
total = mn*xn;
dt = 1e-3; %time step
t = 0:dt:40;
m = linspace(0,1,mn);
dm = m(2) - m(1);
x = linspace(0,40,xn);
dx = x(2) - x(1);
[X,M] = meshgrid(x,m);
tn = length(t);
m_fine = [linspace(0,0.1,100) linspace(0.1,0.9,100) linspace(0.9,1,100)];



%define activation modulus, signal factor, and response to signal factor
% s = @(t) (1+sin(t));

% s = @(t) 1+sin(t);

alpha = 0.5;
beta = 1;
gamma = .5;

f = @(s) beta*(s-1);

[g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example3(alpha,beta,gamma);
IC_1_d_m = IC_uniform(.05,.35);

% g = @(m) alpha*m.*(1-m);%(1-m)/4;
% sigma_inv = @(t,m) (m.*exp(alpha.*t))./((1-m)+m.*exp(alpha*t));
% sigma = @(m2,m1) log((m2.*(1-m1))./((1-m2).*(m1)))/alpha;
% IC_1_d_m = @(m) 10/3*(m>=0.05).*(m<=0.35);

Soln = @(t,s) g(sigma_inv(-int_f_s(t),s))./(g(s)).*IC_1_d_m(sigma_inv(-int_f_s(t),s));


%find x locations where D large
D_cut = .5;
lambda_cut = 0.5;

%parameter values
D_large = 1;
D_small = D_large*1e-2;

lambda_large = 0.25;
lambda_small = lambda_large*1e-2;

%nonautonomous diffusion, proliferation


D_nonaut = @(t) D_large + (D_small - D_large)*uniform_cdf(0.05,0.35,psi(t));
lambda_nonaut = @(t) lambda_small + (lambda_large - lambda_small)*uniform_cdf(0.05,0.35,psi(t));

% 
% 
%initial condition
% u0 = 10/3;
% IC = u0*(X<=6).*(X>=2).*(M>=.2).*(M<=.4);
IC = (X<5).*IC_1_d_m(M);
IC1d = sum(IC)/mn;
IC = IC(:);
% 

%initialize
u = zeros(total,tn);
u(:,1) = IC;

tic

%and nonautonomous
z_nonaut = RD_sim_nonaut_ex1(D_nonaut,lambda_nonaut,t,x,IC1d);

toc

umax = max(max(u));


%     figure%('units','normalized','outerposition',[0 0 1 1]);
%     hold on
% 
%     last_period = find(t>t(end)-6*pi/gamma,1,'first');
% 
%     %silly trick for making legend -- plot on top of self.
%       plot(x,-ones(length(x),1),'b--','linewidth',.5)
%     plot(x,-ones(length(x),1),'r-','linewidth',.5)
% 
% 
%     
%     legend('Proliferating','Diffusing')
% 
%     
%     
%     for i = last_period:1500:tn
%         if psi(t(i)) > .20
%             m_col = 'b--';
%         else
%             m_col = 'r-';
%         end
% 
%         plot(x,z_nonaut(i,:),m_col,'linewidth',.5)
%     end
% 
%     xlabel('Location (x)')
%     ylabel('w(t,x)')
% 
% %     set(gca,'xtick',[])
%     
%     axis([5 40 0 1.1])
%     
%     arrow([12 .5],[30 .5])
%     text(25,.6,'time increasing','fontsize',30)
% 
%     title('Averaged Fishers Equation, Example 3','interpreter','latex')
% 
%     text(1,1.16,'(a)')
%     
%     scale = 0.05;
% 
%     pos = get(gca,'position');
%     pos(2) = pos(2) +scale*pos(4);
%     pos(4) = (1-scale)*pos(4);
%     set(gca,'position',pos)
%     
%     
%     
%     exportfig(gcf,['nonaut_fisher_profile_ex3.eps'])
%     saveas(gcf,['nonaut_fisher_profile_ex3.fig'])


figure
hold on

contourf(x,t(1:100:tn),z_nonaut(1:100:tn,:),20,'edgecolor','none')

% %find where switch occurs for more than half population diffusing
switchtime = find(abs(uniform_cdf(.05,.35,psi(t))-.5)<1e-3);
for i = 1:length(switchtime)
    plot([-5 30],[t(switchtime(i)) t(switchtime(i))],'k')
end

% reg1 = sigma(.5,.05);
% reg2 = sigma(.5,.35);

% plot([-5 30],[t(switchtime(1)) t(switchtime(1))],'k')
% plot([-5 30],[t(switchtime(2)) t(switchtime(2))],'k')
% plot([-5 30],[reg1 reg1],'k')
% plot([-5 30],[reg2 reg2],'k')

xlabel('Location ($x$)','interpreter','latex')
ylabel('Time ($t$)','interpreter','latex')

%make contour blue
map1 = ones(500,1); 
unitgrad = linspace(0,1,500);
map = [flipud(unitgrad') flipud(unitgrad') map1 ];
colormap(map)
caxis([0 1])
colorbar

axis([0 25 0 40])

title('$w(t,x)$ for Example 3','interpreter','latex')

text(-4,42,'(a)','fontsize',15)

% text(16,22.5,'Primarily proliferating')
% text(16,7.5,'Primarily diffusing')
% text(16,.75,'Primarily proliferating')

period = 2*pi/gamma;

text(25,2 ,'P','horizontalalignment','right')
text(25,6,'D','horizontalalignment','right')

for i = 1:4
    text(25,12 + (i-1)*period,'P','horizontalalignment','right')
    text(25,19 + (i-1)*period,'D','horizontalalignment','right')
end

exportfig(gcf,'nonaut_fisher_profile_ex3.eps','color','rgb','renderer','opengl','fontsize',1.2)
saveas(gcf,'nonaut_fisher_profile_ex3.fig')





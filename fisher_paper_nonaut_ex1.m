%fisher_paper_nonaut_ex2.m written 11-14-16 by JTN to simulate
%nonautonomous Fisher's Equation for example 2

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


alpha = 0.1;
[g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example1(alpha);


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
IC1d = double(x<5);
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


    figure%('units','normalized','outerposition',[0 0 1 1]);
    hold on

 %silly trick for making legend -- plot OFF AXIS.
    plot(x,-ones(length(x),1),'b--','linewidth',.5)
    plot(x,-ones(length(x),1),'r-','linewidth',.5)

    legend('Proliferating','Diffusing')
    
    for i = 1:6000:tn
        if psi(t(i)) > .20
            m_col = 'b--';
        else
            m_col = 'r-';
        end

        plot(x,z_nonaut(i,:),m_col,'linewidth',.5)
    end

    xlabel('Location (x)')
    ylabel('w(t,x)')

    axis([0 30 0 1.05])
    arrow([4.5 .15],[24,.15])
    text(18,.2,'time increasing','fontsize',30)
    
    text(-4,1.1,'(a)','fontsize',25)
    title('Averaged Fishers Equation, Example 1','interpreter','latex')


    scale = 0.05;

    pos = get(gca,'position');
    pos(2) = pos(2) +scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca,'position',pos)
    
    exportfig(gcf,['nonaut_fisher_profile_ex1.eps'])
    saveas(gcf,['nonaut_fisher_profile_ex1.fig'])


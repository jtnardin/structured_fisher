%sim_nonaut_w_tc_range_ex1.m written 11-15-16 by JTN to simulate
%nonautonomous fisher equation over a range of parameter values

clear all; clc

%define mesh
xn = 301; %number of x points
dt = 1e-3; %time step
t = 0:dt:40;
x = linspace(0,40,xn);

alpha = linspace(0,.2,50);

z_nonaut = cell(length(alpha),1);
IC1d = double(x<=5);


for i = 1:length(alpha)
   
    %load functions
    [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example1(alpha(i));
    
    %parameter values
    D_large = 1;
    D_small = D_large*1e-2;

    lambda_large = 0.25;
    lambda_small = lambda_large*1e-2;

    %nonautonomous diffusion, proliferation rates

    D_nonaut = @(t) D_large + (D_small - D_large)*uniform_cdf(0.05,0.35,psi(t));
    lambda_nonaut = @(t) lambda_small + (lambda_large - lambda_small)*uniform_cdf(0.05,0.35,psi(t));
   
    
    %sim nonaut version
    nonaut_soln = RD_sim_nonaut_ex1(D_nonaut,lambda_nonaut,t,x,IC1d);
    
    
%     just get final profile
    z_nonaut{i} = nonaut_soln(end,:);
end




%%%%make a contour



z_cont = zeros(length(alpha),length(x));

for i = 1:length(alpha)
    z_cont(i,:) = z_nonaut{i};
end

figure
hold on

contourf(x,alpha,z_cont,20,'edgecolor','none')


%plot different regions
reg2 = log(13/7)/40;
reg1 = log(19)/40;

plot([-5 30],[reg2 reg2],'k')
plot([-5 30],[reg1 reg1],'k')

axis([0 25 alpha(1) alpha(end) 0 1])
view(2)
xlabel('Location ($x$)','fontsize',30,'interpreter','latex')
ylabel('Growth Rate ($\alpha$)','fontsize',30,'interpreter','latex')

%make contour blue
map1 = linspace(1,.5,500)'; 
unitgrad = linspace(0,1,500);
map = [flipud(unitgrad')  map1  flipud(unitgrad')];
colormap(map)
caxis([0 1])
colorbar

title('$w(t=40,x)$ for various values of $\alpha$','fontsize',30,'interpreter','latex')

text(-4.7,.21,'(b)','fontsize',15)
text(25,.135,'Entire Activation','horizontalalignment','right')
text(25,.04,'Activation','horizontalalignment','right')
text(25,0.0075,'No Activation','horizontalalignment','right')


exportfig(gcf,'Ex1_nonaut_range_alpha.eps','color','rgb','resolution',600,'renderer','opengl','fontsize',1.2)
saveas(gcf,'Ex1_nonaut_range_alpha.fig')


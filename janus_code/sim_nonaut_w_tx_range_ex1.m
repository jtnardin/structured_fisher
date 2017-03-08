%sim_nonaut_w_tc_range_ex1.m written 11-15-16 by JTN to simulate
%nonautonomous fisher equation over a range of parameter values

clear all; clc

%define mesh
xn = 301; %number of x points
dt = 1e-3; %time step
t = 0:dt:40;
x = linspace(0,40,xn);

alpha = [0.01 0.1:.1:.9];

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

fig_final = figure('units','normalized','outerposition',[0 0 1 1]);
hold on

% c = distinguishable_colors(length(alpha));

for i = 1:length(alpha)
    
% %     plot(x,z{i}(:,end),'color',c(i,:),'linewidth',6)
%     if i == 2
%       plot(x,z_nonaut{i}(:),'--','color',c(i,:),'linewidth',6)
% 
%     else
      plot(x,z_nonaut{i}(:),'linewidth',6)
%     end
end
% 
% matrix_legend = cell(param_length,1);
% 
% for i = 1:param_length
%    matrix_legend{i} = num2str(param(i)); 
% end
% 
% h=legend(matrix_legend);
% set(h,'fontsize',30)
axis([0 20 0 1.05])

text(-1.5,1.1,'(d)','fontsize',30)

xlabel('x','fontsize',35)
ylabel(['w(t=' num2str(t(end)) ',x)'],'fontsize',35)

title('Example 2 for various $\beta$ values','interpreter','latex','fontsize',35)

set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])

export_fig(gcf,'Ex2_nonaut_range_beta_transition.eps')
saveas(gcf,'Ex2_nonaut_range_beta_transition.fig')
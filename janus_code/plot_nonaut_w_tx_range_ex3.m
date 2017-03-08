clear all; clc

load('ex3_range_gamma_transition_11_7')

param = gamma; %what parameter was varied
param_length = length(gamma)+1; %floor(length(beta)/2)

dgamma = gamma(3) - gamma(2);

param(1) = gamma(2) - dgamma;
param(end) = gamma(end-1)+dgamma;
param(end+1) = gamma(end) + dgamma;

z_nonaut = cell(param_length,1);

IC1d = double(x<=5);

for i = 1:param_length
   
    [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example3(alpha,beta,param(i));
    
     %parameter values
    D_large = 1;
    D_small = D_large*1e-2;

    lambda_large = 0.25;
    lambda_small = lambda_large*1e-2;

    %nonautonomous diffusion, proliferation

    D_nonaut = @(t) D_large + (D_small - D_large)*uniform_cdf(0.05,0.35,psi(t));
    lambda_nonaut = @(t) lambda_small + (lambda_large - lambda_small)*uniform_cdf(0.05,0.35,psi(t));
   
    nonaut_soln = RD_sim_nonaut_ex1(D_nonaut,lambda_nonaut,t,x,IC1d);
    
    z_nonaut{i} = nonaut_soln(end,:);
end

figure
hold on

plot(x,-ones(length(x),1),'k-')
plot(x,-ones(length(x),1),'b--')
plot(x,-ones(length(x),1),'r-.')

h= legend('No activation','Some activation','Full activation','location','southwest');
% set(h,'fontsize',30)


for i = 1:param_length
    
    
  if param(i) >= 1.615
      color_m = 'k-';
  elseif param(i) > 0.34
      color_m = 'b--';
  else
      color_m = 'r-.';
  end
    
  plot(x,z_nonaut{i}(:),color_m,'linewidth',1)
    
    
end

axis([0 20 0 1.05])

text(-1.5,1.1,'(d)')%,'fontsize',30)

xlabel('Location (x)')%,'fontsize',35)
ylabel(['w(t=' num2str(t(end)) ',x)'])%,'fontsize',35)

title('Example 3 for various $\gamma$ values','interpreter','latex')%,'fontsize',35)

set(gca,'fontsize')%,30)
% set(gcf,'color',[1 1 1])

exportfig(gcf,'Ex3_nonaut_range_beta_transition.eps')
saveas(gcf,'Ex3_nonaut_range_beta_transition.fig')
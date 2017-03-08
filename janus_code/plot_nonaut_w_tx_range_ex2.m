clear all; clc

load('ex2_range_beta_transition_11_7')


param = beta; %what parameter was varied

param(1) = 2.1029;
param(end) = 6.1271;

%add more terms past entire activationthreshold
param(end+1) = 2*param(end) - param(end-1);
param(end+1) = 2*param(end) - param(end-1);
param(end+1) = 2*param(end) - param(end-1);

param_length = length(param); %floor(length(beta)/2)

z_nonaut = cell(param_length,1);

IC1d = double(x<=5);

for i = 1:param_length
   
    [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example2(alpha,param(i),gamma);
    
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

figure%('units','normalized','outerposition',[0 0 1 1]);
hold on

plot(x,-ones(length(x),1),'k-')
plot(x,-ones(length(x),1),'b--')
plot(x,-ones(length(x),1),'r-.')

h= legend('No activation','Some activation','Full activation','location','southwest');
% set(h,'fontsize',30)


for i = 1:param_length
  
  if param(i) <= 2.55
      color_m = 'k-';
  elseif param(i) <= 5.68
      color_m = 'b--';
  else
      color_m = 'r-.';
  end
    
  plot(x,z_nonaut{i}(:),color_m,'linewidth',1)
    
end

axis([0 20 0 1.05])

text(-1.5,1.1,'(b)')%,'fontsize',30)

xlabel('Location (x)')%,'fontsize',30)
ylabel(['w(t=' num2str(t(end)) ',x)'])%,'fontsize',35)

title('Example 2 for various $\beta$ values','interpreter','latex')%,'fontsize',35)
arrow([5 .5],[17 .5],'color','k')
text(2,.55,'$\beta$ increasing','interpreter','latex')


% set(gca,'fontsize',30)
% set(gcf,'color',[1 1 1])


scale = 0.05;

pos = get(gca,'position');
pos(2) = pos(2) +scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca,'position',pos)

exportfig(gcf,'Ex2_nonaut_range_beta_transition.eps')
saveas(gcf,'Ex2_nonaut_range_beta_transition.fig')

clear all; clc

% load('ex3_range_gamma_transition_11_7')
% 
% param = gamma; %what parameter was varied
% param_length = length(gamma)+1; %floor(length(beta)/2)
% 
% dgamma = gamma(3) - gamma(2);
% 
% param(1) = gamma(2) - dgamma;
% param(end) = gamma(end-1)+dgamma;
% param(end+1) = gamma(end) + dgamma;


%define mesh
xn = 301; %number of x points
dt = 1e-3; %time step
t = 0:dt:40;
x = linspace(0,40,xn);

alpha = 0.5;
beta = 1;
gamma = [linspace(0,.1,10) linspace(.101,.9,50) linspace(.91,1.9,20)];

param_length = length(gamma);

z_nonaut = cell(param_length,1);

IC1d = double(x<=5);

for i = 1:param_length
   
    if gamma(i) == 0 %special case for gamma = 0. 
    
        psi = @(t) .5;
        
    else
    
        [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example3(alpha,beta,gamma(i));
    
    end
        
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

% figure
% hold on
% 
% plot(x,-ones(length(x),1),'k-')
% plot(x,-ones(length(x),1),'b--')
% plot(x,-ones(length(x),1),'r-.')
% 
% h= legend('No activation','Some activation','Full activation','location','southwest');
% % set(h,'fontsize',30)
% 
% 
% for i = 1:param_length
%     
%     
%   if param(i) >= 1.615
%       color_m = 'k-';
%   elseif param(i) > 0.34
%       color_m = 'b--';
%   else
%       color_m = 'r-.';
%   end
%     
%   plot(x,z_nonaut{i}(:),color_m,'linewidth',1)
%     
%     
% end
% 
% axis([0 20 0 1.05])
% 
% text(-1.5,1.1,'(d)')%,'fontsize',30)
% 
% xlabel('Location (x)')%,'fontsize',35)
% ylabel(['w(t=' num2str(t(end)) ',x)'])%,'fontsize',35)
% 
% title('Example 3 for various $\gamma$ values','interpreter','latex')%,'fontsize',35)
% 
% set(gca,'fontsize')%,30)
% % set(gcf,'color',[1 1 1])





z_cont = zeros(length(gamma),length(x));

for i = 1:length(gamma)
    z_cont(i,:) = z_nonaut{i};
end
gamma_crit1 = .34;
gamma_crit2 = 1.615;

figure
hold on

contourf(x,gamma,z_cont,20,'edgecolor','none')

plot([-5 30],[gamma_crit1 gamma_crit1],'k')
plot([-5 30],[gamma_crit2 gamma_crit2],'k')

axis([0 25 gamma(1) gamma(end)])
view(2)
xlabel('Location ($x$)','interpreter','latex')
ylabel('Treatment Frequency ($\gamma$)','interpreter','latex')

%make contour blue
map1 = linspace(1,.5,500)'; 
unitgrad = linspace(0,1,500);
map = [flipud(unitgrad') map1 flipud(unitgrad')];
colormap(map)
caxis([0 1])
colorbar

title('$w(t=40,x)$ for various values of $\gamma$','interpreter','latex')

text(-4.2,2,'(b)','fontsize',15)

text(25,1,'Activation','horizontalalignment','right')
text(25,.08,'Entire Activation','horizontalalignment','right')
text(25,1.75,'No Activation','horizontalalignment','right')


exportfig(gcf,'Ex3_nonaut_range_beta_transition.eps','color','rgb','renderer','opengl','fontsize',1.3)
saveas(gcf,'Ex3_nonaut_range_beta_transition.fig')


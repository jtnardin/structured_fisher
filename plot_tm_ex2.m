clear all; clc

t = linspace(0,15,300);
s = linspace(0,1,400);

[T,S] = meshgrid(t,s);
alpha = 0.5;

[g,sigma,sigma_inv,s_t,f,int_f_s] = g_sigma_h_example2(alpha);

phi = IC_uniform(.05,.35);

Soln = @(t,s) g(sigma_inv(-int_f_s(t),s))./(g(s)).*phi(sigma_inv(-int_f_s(t),s));

A = Soln(T,S);
A(S==0)=phi(s(1));
A(S==1)=phi(s(end));

B = log(A);
B(isinf(B))=NaN;

figure
hold on

contourf(T,S,B,'edgecolor','none')

plot(t,sigma_inv(int_f_s(t),0.35),'color','k','linewidth',1)
plot(t,sigma_inv(int_f_s(t),0.15),'color','k','linewidth',1)
plot(t,sigma_inv(int_f_s(t),0.05),'color','k','linewidth',1)


%plot m = 0.5

plot([-2 17],[.5 .5],'k--','linewidth',0.5)
text(1.7,.3,'Inactive Population')
text(1.7,.6,'Active Population')

axis([0 15 0 1])

title('Activation Profile for Example 2','interpreter','latex')
xlabel(' Time (t) ')
ylabel(' Activation level (m) ')


view(2)


map1 = ones(500,1); 
unitgrad = linspace(0,1,500);

map = [ map1 flipud(unitgrad') flipud(unitgrad') ];


colormap(map)

caxis([min(min(log(A(A~=0)))) max(max(log(A)))])
colorbar

exportfig(gcf,'ex2_tm.eps','color','rgb','renderer','opengl','fontsize',1.2)
saveas(gcf,'ex2_tm.fig')

t = linspace(0,5,300);
s = linspace(0,1,400);

[T,S] = meshgrid(t,s);

alpha = 1;

[g,sigma,sigma_inv,s_t,f,int_f_s] = g_sigma_h_example1(alpha);

phi = IC_uniform(0.05,0.35);

Soln = @(t,s) g(sigma_inv(-t,s))./(g(s)).*phi(sigma_inv(-t,s));

A = Soln(T,S);
%make white in plot
A(A==0)=NaN;

A(S==0)=phi(s(1));
A(S==1)=phi(s(end));

figure(1)
hold on

B = log(A);

B(isinf(B))=NaN;

contourf(T,S,B,'edgecolor','none')

%plot characteristics
plot(t,sigma_inv(int_f_s(t),0.35),'color','k','linewidth',1)
plot(t,sigma_inv(int_f_s(t),0.15),'color','k','linewidth',1)
plot(t,sigma_inv(int_f_s(t),0.05),'color','k','linewidth',1)


%plot m = 0.5

plot([-2 7],[.5 .5],'k--','linewidth',0.5)
text(5.1,.45,'Inactive Population','horizontalalignment','right')
text(5.1,.55,'Active Population','horizontalalignment','right')

title('Activation Profile for Example 1','interpreter','latex')
xlabel('Time (t) ')
ylabel('Activation level (m)')

axis([0 5 0 1])

caxis([min(min(log(A(A~=0)))) max(max(log(A)))]/1.5)

map1 = ones(500,1); 
unitgrad = linspace(0,1,500);

map = [ map1 flipud(unitgrad') flipud(unitgrad') ];


colormap(map)
colorbar
view(2)

exportfig(gcf,'ex1_tm.eps','color','rgb','renderer','opengl','fontsize',1.2)
saveas(gcf,'ex1_tm.fig')

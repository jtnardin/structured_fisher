t = linspace(0,13,300);
s = linspace(0,1,400);

[T,S] = meshgrid(t,s);

alpha = .5;

[g,sigma,sigma_inv,s_t,f,int_f_s] = g_sigma_h_example3(alpha);

phi = IC_uniform(0.05,0.35);

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
text(3,0,'Subpopulation 1','fontsize',20)
text(3,1.1,'Subpopulation 2','fontsize',20)

title('Activation Distribution over time for Example 2','interpreter','latex')
xlabel(' Time (t) ')
ylabel('Activation level (m) ')
yticks([0 0.2 0.4 0.6 0.8 1])


axis([0 12 -.1 1.2])

map1 = ones(500,1); 
unitgrad = linspace(0,1,500);

map = [flipud(unitgrad') flipud(unitgrad') map1 ];

colormap(map)
colorbar
caxis([min(min(log(A(A~=0)))) max(max(log(A)))])

view(2)

exportfig(gcf,'ex3_tm.eps')
saveas(gcf,'ex3_tm.fig')

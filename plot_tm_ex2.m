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


figure
hold on

contourf(T,S,B,'edgecolor','none')

plot(t,sigma_inv(int_f_s(t),0.35),'color',[0.5 0.5 0.5],'linewidth',1)

% h=legend('Distribution of m','$h(t;\underline{m})$');
% set(h,'interpreter','latex')

plot(t,sigma_inv(int_f_s(t),0.15),'color',[0.5 0.5 0.5],'linewidth',1)
plot(t,sigma_inv(int_f_s(t),0.05),'color',[0.5 0.5 0.5],'linewidth',1)

% caxis([min(min(log(A(A~=0)))) max(max(log(A)))])


title('Distribution along $m$ over time in Example 2','interpreter','latex')
xlabel('t')
ylabel('m')


view(2)

caxis([min(min(log(A(A~=0)))) max(max(log(A)))])
colorbar

exportfig(gcf,'ex2_tm.eps')
saveas(gcf,'ex2_tm.fig')

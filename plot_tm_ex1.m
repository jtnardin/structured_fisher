t = linspace(0,5,300);
s = linspace(0,1,400);

[T,S] = meshgrid(t,s);

alpha = 1/2;

g = @(s) alpha*s.*(1-s);

sigma_inv = @(t,m) (m.*exp(alpha*t))./((1-m)+m.*exp(alpha*t));

phi = @(m) 10/3*(m>=0).*(m<=0.3);

Soln = @(t,s) g(sigma_inv(-t,s))./(g(s)).*phi(sigma_inv(-t,s));

A = Soln(T,S);

figure

contourf(T,S,A,'edgecolor','none')

title('$\sigma^{-1}(t;\underline{m})$, Example 1','interpreter','latex')
xlabel('t')
ylabel('m')

colorbar
view(2)

exportfig(gcf,'ex1_tm.eps')
saveas(gcf,'ex1_tm.fig')

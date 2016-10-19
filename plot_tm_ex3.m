t = linspace(0,15,300);
s = linspace(0,1,400);

[T,S] = meshgrid(t,s);

alpha = .5;

g = @(s) alpha*s.*(1-s);

sigma_inv = @(t,m) (m.*exp(alpha*t))./((1-m)+m.*exp(alpha*t));

int_f_s = @(t) 2*(1-cos(t));

phi = @(m) 10/3*(m>=0).*(m<=0.3);

Soln = @(t,s) g(sigma_inv(-int_f_s(t),s))./(g(s)).*phi(sigma_inv(-int_f_s(t),s));

A = Soln(T,S);

figure

contourf(T,S,A,'edgecolor','none')

title(' Example 3 $\sigma^{-1}(t;\underline{m})$','interpreter','latex')
xlabel('t')
ylabel('m')

colorbar
view(2)

exportfig(gcf,'ex3_tm.eps')
saveas(gcf,'ex3_tm.fig')

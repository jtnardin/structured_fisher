t = linspace(0,1,300);
s = linspace(1,3,80);

[T,S] = meshgrid(t,s);

g = @(s) sqrt(s-1).*(3-s);

sigma_inv = @(t,s) 1+2*(tanh(.5*(sqrt(2)*t + 2*atanh(sqrt(s-1)/sqrt(2))))).^2;

phi = @(s) 10*(s-1.1).^2.*(1.4-s).^2.*(s>=1.1).*(s<=1.4) + (s-2).^2.*(2.5-s).^2.*(s>=2).*(s<=2.5);

Soln = @(t,s) g(sigma_inv(-t,s))./(g(s)).*phi(sigma_inv(-t,s)).*(sigma_inv(t,1)<s);

A = Soln(T,S);

figure

surf(T,S,A,'edgecolor','none')

xlabel('t')
ylabel('s')

caxis([0,0.01])
colorbar

axis([0 1 1 3 0 0.013])

view(2)
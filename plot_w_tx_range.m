% clear all; clc
% 
% load('ex1_range_alpha_2')

x = linspace(0,60,2*xn-1);

figure
hold on

for i = 1:length(alpha)
    
    plot(x,z{i}(:,end),'linewidth',.25)
    
end


arrow([10,0.4],[30,0.4])
text(10,0.45,'$\alpha$ increasing','interpreter','latex')

axis([0 40 0 1])

xlabel('x')
ylabel('w(t,x)')

title('Example 1, w(t=30,x) for various $\alpha$ values','interpreter','latex')

exportfig(gcf,'Ex1_range_alpha.eps')
saveas(gcf,'Ex1_range_alpha.fig')
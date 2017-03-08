clear all; clc

load('ex2_range_beta_transition_11_7')

fig_final = figure('units','normalized','outerposition',[0 0 1 1]);
hold on

param = beta; %what parameter was varied
param_length = length(beta); %floor(length(beta)/2)


c = distinguishable_colors(param_length);

for i = 1:param_length
    
    plot(x,z{i}(:,end),'color',c(i,:),'linewidth',6)
    
end

matrix_legend = cell(param_length,1);

for i = 1:param_length
   matrix_legend{i} = num2str(param(i)); 
end

h=legend(matrix_legend);
set(h,'fontsize',30)
axis([0 20 0 1.05])

text(-1.5,1.1,'(d)','fontsize',30)

xlabel('x','fontsize',35)
ylabel(['w(t=' num2str(t(end)) ',x)'],'fontsize',35)

title('Example 2 for various $\beta$ values','interpreter','latex','fontsize',35)

set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])

export_fig(gcf,'Ex2_range_beta_transition.eps')
saveas(gcf,'Ex2_range_beta_transition.fig')
clear all; clc

load('ex3_range_gamma_transition_11_7')

% c = distinguishable_colors(floor(length(beta)/2));

param = gamma;
param_length = length(gamma);

colors = distinguishable_colors(param_length);

step=floor(tn/4);
count = 1;
for j = step:step:tn

    fig_final = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    for i = 1:length(param)
                
        plot(x,z{i}(:,j),'color',colors(i,:),'linewidth',6)
    
    end

    title([' t = ' num2str(t(j))])
    count = count + 1;
    
end

text(-2,1.1,'(b)','fontsize',35)

matrix_legend = cell(length(param),1);

for i = 1:length(param)
   matrix_legend{i} = num2str(param(i)); 
end

h=legend(matrix_legend);
set(h,'fontsize',30)
axis([0 25 0 1.1])

set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])

xlabel('x','fontsize',35)
ylabel('w(t=30,x)','fontsize',35)

title('Example 3 for various $\gamma$ values','interpreter','latex','fontsize',35)

export_fig(gcf,'Ex3_range_gamma.eps')
saveas(gcf,'Ex3_range_gamma.fig')



figure
hold on

[g,sigma,sigma_inv,s,f,int_f_s] = g_sigma_h_example3(alpha,beta,gamma(1));

h = @(t) sigma_inv(-int_f_s(t),1/2);


last_period = find(t>t(end)-4*pi,1,'first');

%silly trick for making legend -- plot on top of self.
plot(x,z{end}(:,last_period),'b','linewidth',1)
plot(x,z{end}(:,last_period),'r','linewidth',1)

for i = last_period:750:tn
    if h(t(i)) > .20
        m_col = 'b';
    else
        m_col = 'r';
    end

    plot(x,z{1}(:,i),m_col,'linewidth',.25)
end

xlabel('x')
ylabel('w(t,x)')

axis([0 20 0 1.1])

title('Example 3, w(t,x)')

legend('Proliferating','Diffusing')

% exportfig(gcf,'nonaut_fisher_profile_ex3.eps')
% saveas(gcf,'nonaut_fisher_profile_ex3.fig')

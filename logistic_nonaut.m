%func logistic_nonaut.m written 5-30-16 by JTN to simulate a nonautonomous
%logistic equation to see if it converges to the logistic equation over
%time

clear all; clc

t = linspace(0,50,100);


%exact autonomous solution
v_star = @(t,v0) v0*exp(t)./(v0*exp(t)+1-v0);

m0 = [0];

a = [1e-4,0.1,.5,.9];

for i = 1:length(a)


    %initial value
    v0 = a(i);

    %simulate nonautonomous ODE
    [t,v] = ode45(@(t,v) nonaut_log_ode(t,v,m0),t,v0);

    figure

    axes('Position',[.005 .005 .99 .99],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
    
    subplot(2,1,1)

    hold on
    
    plot(t,v_star(t,v0),'b-')
    plot(t,v,'b--')

    xlabel('Time')
    ylabel('v(t)')
    title(['autonomous vs. nonautonomous solutions, $\underbar{m}$ = ' num2str(m0) ' v0 = ' num2str(v0)],'interpreter','latex')

    if i == 1
       legend('v^*(t)','v(t)','location','southeast')         
    end
    
    axis([0 50 0 1.6])
    
    subplot(2,1,2)

    plot(t,v_star(t,v0)-v,'b.')


    xlabel('Time')
    ylabel('v^*(t) - v(t)')
    title('Difference between solutions')

    axis([0 50 -0.8 0])
    
    set(gcf,'color',[1 1 1])
    
    
    export_fig(gcf,['log_nonaut_m_' num2str(m0) '_' num2str(i) '.eps'])
    
end
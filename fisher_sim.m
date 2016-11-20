t = linspace(0,30,1000);
x = linspace(-10,10,100);
u = RD_sim(.1,.25,t,x,double(x<=-5));



for i = 1:6
    hold on
    plot(x,u(400+100*i,:))
end

xlabel('x')
ylabel('u(t,x)')

title('Numerical simulation of Fishers Equation')



axis([-7 8 0 1.1])

arrow([-4,.45],[4,.45])
text(2,.5,'t increasing','fontsize',30)

% exportfig(gcf,'fisher_sim.eps')

figure
hold on

plot(x,double(x<=-5));

axis([-10 10 0 1.1])

text(0,0.4,'Wound space','fontsize',30)
text(-9.5,1.05,'wound margin','fontsize',30)


xlabel('x')
ylabel('u(0,x)')
title('Initial condition')

exportfig(gcf,'IC_fisher.eps')

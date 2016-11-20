clear all; clc

alpha = .1:.1:1;

z = cell(length(alpha),1);


for i = 1:length(alpha)
    i
    z{i} = fisher_paper_ex1_f(alpha(i));
end

mn = 81; %number of m points
xn = 151; %number of x points
total = mn*xn;
dt = 1e-3; %time step
t = 0:dt:15;
m = linspace(0,1,mn);
dm = m(2) - m(1);
x = linspace(0,30,xn);
dx = x(2) - x(1);
[X,M] = meshgrid(x,m);
tn = length(t);
m_fine = [linspace(0,0.1,100) linspace(0.1,0.9,100) linspace(0.9,1,100)];

save('ex1_range_alpha_2.mat');
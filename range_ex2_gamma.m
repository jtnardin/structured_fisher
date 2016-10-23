clear all; clc

alpha = 0.5;
beta = 3;

gamma_inv = 1:2:8;

gamma = -1./gamma_inv;

y = cell(length(gamma),1);
z = cell(length(gamma),1);


for i = 1:length(gamma)
    [y{i},z{i}] = fisher_paper_ex2_f(alpha,beta,gamma(i));
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

save('ex2_range_gamma.mat');
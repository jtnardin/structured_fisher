clear all; clc

alpha = .1:.1:1;

z = cell(length(alpha),1);


% Just input time, space and m grids, IC?  


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


IC_1_d_m = IC_uniform(.05,.35);%@(m) 1*(m>=0).*(m<=0.3);

%initial condition
IC = IC_1_d_m(M).*(X<=5);%.*exp(-M);
IC1d = sum(IC)/mn;
IC = IC(:);

for i = 1:length(alpha)
    i
    z{i} = fisher_paper_ex1_f_jan(alpha(i),t,x,m,X,M,IC);
end


save('/lustre/janus_scratch/jona8898/fisher_struct/ex1_range_alpha.mat');
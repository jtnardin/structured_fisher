%%%g_sigma_h_example1.m written 10-19-16 by JTN to give all functions
%%%necessary for example 1 in structured fisher paper
function [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example1(alpha)

%     alpha = 0.5;
    
    g = @(m) alpha*m.*(1-m);%(1-m)/4;
    sigma_inv = @(t,m) (m.*exp(alpha*t))./((1-m)+m.*exp(alpha*t));
    sigma = @(m2,m1) log((m2.*(1-m1))./((1-m2).*(m1)))/alpha;


    s = @(t) 1+sin(t);

    f = @(s) 1;
    
    int_f_s = @(t) t;
    
    psi = @(t) 1./(1+exp(alpha*t));
end
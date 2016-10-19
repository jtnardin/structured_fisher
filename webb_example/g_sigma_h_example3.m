%%%g_sigma_h_example3.m written 10-19-16 by JTN to give all functions
%%%necessary for example 3 in structured fisher paper
function [g,sigma,sigma_inv,s,f,int_f_s] = g_sigma_h_example3(alpha)


    g = @(m) alpha*m.*(1-m);

    sigma = @(m2,m1) log((m2.*(1-m1))./((1-m2).*(m1)))/alpha;
    sigma_inv = @(t,m) (m.*exp(alpha*t))./((1-m)+m.*exp(alpha*t));


    beta = 4;

    
    f = @(s) beta*(s-1);
    s = @(t) 1 + sin(t);
    
    int_f_s = @(t) beta*(1-cos(t));
    
end
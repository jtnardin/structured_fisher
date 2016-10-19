%%%g_sigma_h_Webb_example.m written 10-19-16 by JTN to give all functions
%%%necessary for example 1.1 in Webb 2008
function [g,sigma,sigma_inv,s,f,int_f_s] = g_sigma_h_WebbExample

    %growth rate
    g = @(s) sqrt(s-1).*(3-s);
    sigma = @(m2,m1) sqrt(2)*atanh(sqrt(m2-1)/sqrt(2)) - sqrt(2)*atanh(sqrt(m1-1)/sqrt(2));
    sigma_inv = @(t,s) 1+2*(tanh(.5*(sqrt(2)*t + 2*atanh(sqrt(s-1)/sqrt(2))))).^2;


    s = @(t) 1+sin(t);

    f = @(s) 1;
    
    int_f_s = @(t) t;
end
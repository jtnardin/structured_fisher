%%%g_sigma_h_example3.m written 10-19-16 by JTN to give all functions
%%%necessary for example 3 in structured fisher paper
function [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example3(alpha,beta,gamma)


    if nargin == 1
        beta=4;
        gamma = 1;
    elseif nargin == 2
        gamma = 1;
    end

    g = @(m) alpha*m.*(1-m);

    sigma = @(m2,m1) log((m2.*(1-m1))./((1-m2).*(m1)))/alpha;
    sigma_inv = @(t,m) (m.*exp(alpha*t))./((1-m)+m.*exp(alpha*t));

    
    f = @(s) beta*(s-1);
    s = @(t) 1 + sin(gamma*t);
    
    int_f_s = @(t) beta/gamma*(1-cos(gamma*t));
    
    psi = @(t) 1./(1+exp(alpha*beta/gamma*(1-cos(gamma*t))));
    
end
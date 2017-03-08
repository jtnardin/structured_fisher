%%%g_sigma_h_example2.m written 10-19-16 by JTN to give all functions
%%%necessary for example 2 in structured fisher paper
function [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example2(alpha,beta,gamma)

    if nargin == 1
        beta = 3;
        gamma = -1/4;
    elseif nargin == 2
        gamma = -1/4;
    end


    g = @(s) alpha*s.*(1-s);

    sigma = @(m2,m1) log((m2.*(1-m1))./((1-m2).*(m1)))/alpha;
    sigma_inv = @(t,m) (m.*exp(alpha*t))./((1-m)+m.*exp(alpha*t));

    
    f = @(s) s-1 ;
    s = @(t) beta*exp(gamma*t);
    
    int_f_s = @(t) beta/gamma*(exp(gamma*t)-1)-t;
    
    psi = @(t) 1./(1+exp(-alpha*t+alpha*beta/gamma*(exp(gamma*t)-1)));
    
end
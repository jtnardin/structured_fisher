function out = exp_cdf(a,b,t)

%     if t < a
%         out = 0;
%     elseif t > b
%         out = 1;
%     else
%         out = t/(b-a);
%     end

    gamma = 1/(exp(-a)-exp(-b));
    out = 0.*(t<a) + gamma*(exp(-a)-exp(-t)).*(t>=a).*(t<=b) + 1.*(t>b) ; 


end
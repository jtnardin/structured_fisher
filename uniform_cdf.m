function out = uniform_cdf(a,b,t)

%     if t < a
%         out = 0;
%     elseif t > b
%         out = 1;
%     else
%         out = t/(b-a);
%     end

    out = 0.*(t<a) + (t-a)/(b-a).*(t>a).*(t<b) + 1.*(t>b) ; 


end
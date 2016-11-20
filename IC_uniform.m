function IC = IC_uniform(a,b)

    IC = @(s) 1/(b-a)*(s>=a).*(s<=b);

end
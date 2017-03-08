function IC = IC_two_humps_Webb

    IC = @(s) 10*(s-1.1).^2.*(1.4-s).^2.*(s>=1.1).*(s<=1.4) + (s-2).^2.*(2.5-s).^2.*(s>=2).*(s<=2.5);

end
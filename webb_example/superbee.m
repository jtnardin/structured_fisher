function phi = superbee(r)

%     phi = max[0,min(2r,1),min(r,2)]

    phi = zeros(size(r));
    
    a = min(2*r,1);
    b = min(r,2);
    phi = max(a,b);

    phi(phi<0)=0;


end
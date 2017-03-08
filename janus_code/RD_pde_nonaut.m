function yprime = RD_pde_nonaut(t,y,f,A)

    yprime = A(t)*y + f(y,t);

end
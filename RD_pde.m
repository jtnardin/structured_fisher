function yprime = RD_pde(t,y,f,A)

    yprime = A*y + f(y);

end
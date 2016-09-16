%function nonaut_log_ode.m written 5-30-16 by JTN as the ODE for
%nonautonomous logistic equation.

function vprime = nonaut_log_ode(t,v,m)

    vprime = v/(t+exp(m)) + v.*(1-v);
%     vprime = -exp(t)*v/m + v.*(1-v);

end
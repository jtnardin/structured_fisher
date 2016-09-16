%negative_sensor.m written 9-16-16 by JTN to compute the sensors for flux
%limiters when the velocity is negative.

function [r_e,r_w,r_w_m1] = negative_sensor(u,xm_int,m_bd_1)

    r_e = (u(xm_int(1:end-1)+1) - u(xm_int(1:end-1)+2))./(u(xm_int(1:end-1)) - u(xm_int(1:end-1)+1));    
    r_w = (u(xm_int) - u(xm_int+1))./(u(xm_int-1) - u(xm_int-i));
    r_e = [r_e,-1];
    %boundary m = 1
    r_w_m1 = (u(m_bd_1(2:end-1)-1) - u(m_bd_1(2:end-1)))./(u(m_bd_1(2:end-1)-2) - u(m_bd_1(2:end-1)-1));
    
    %eliminate NaN values (0/0 -- not steep!)
    r_e(isnan(r_e)) = 1;
    r_w(isnan(r_w)) = 1;
    r_w_m1(isnan(r_w_m1)) = 1;
    
    %set inf values to large value
    r_e(isinf(r_e)) = 100;
    r_w(isinf(r_w)) = 100;
    r_w_m1(isinf(r_w_m1)) = 100;
    
end
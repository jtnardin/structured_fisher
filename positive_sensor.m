%positive_sensor.m written 9-16-16 by JTN to compute the sensors for flux
%limiters when the velocity is positive.

function [r_e,r_w,r_w_m1,r_e_m0] = positive_sensor2(u,xm_int,m_bd_1,m_bd_0,m_bd_0_int,m_bd_nm1_int)

    r_e = (u(xm_int) - u(xm_int-1))./(u(xm_int+1) - u(xm_int));    
    r_w = (u(xm_int(2:end)-1) - u(xm_int(2:end)-2))./(u(xm_int(2:end)) - u(xm_int(2:end)-1));
    r_w = [-1;r_w];
    
    %fix points that sampled x_n-2
    
    r_w(m_bd_0_int) = -1;
    
  
%     boundary
    r_w_m1 = r_e(m_bd_nm1_int);
    r_e_m0 = r_w(m_bd_0_int);
    
    
    %eliminate NaN values (0/0 -- not steep!)
    r_e(isnan(r_e)) = 1;
    r_w(isnan(r_w)) = 1;
    r_w_m1(isnan(r_w_m1)) = 1;
    r_e_m0(isnan(r_e_m0)) = 1;
    
    %set inf values to large value
    r_e(isinf(r_e)) = 100;
    r_w(isinf(r_w)) = 100;
    r_w_m1(isinf(r_w_m1)) = 100;
    r_e_m0(isinf(r_e_m0)) = 100;
    
end
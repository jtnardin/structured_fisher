%fisher_paper_ex2_f.m written 10-22-16 by JTN to simulate
%u_t + (f(s(t))g(m)u)_m = D(m)*u_xx + lambda(m)*u(1-w) for a given value of
%alpha

%Example 1 : Threshold at m = 1 : f(s) = 1 , g(m) = m(1-m)

function [y,z] = fisher_paper_ex2_f(alpha,beta,gamma)


    %Construct vectors of independent variables
    mn = 81; %number of m points
    xn = 151; %number of x points
    total = mn*xn;
    dt = 1e-3; %time step
    t = 0:dt:15;
    m = linspace(0,1,mn);
    dm = m(2) - m(1);
    x = linspace(0,30,xn);
    dx = x(2) - x(1);
    [X,M] = meshgrid(x,m);
    tn = length(t);
    m_fine = [linspace(0,0.1,100) linspace(0.1,0.9,100) linspace(0.9,1,100)];

    [g,sigma,sigma_inv,s,f,int_f_s] = g_sigma_h_example2(alpha,beta,gamma);

    IC_1_d_m = IC_uniform(.05,.35);%@(m) 1*(m>=0).*(m<=0.3);
    Soln = @(t,s) g(sigma_inv(-t,s))./(g(s)).*IC_1_d_m(sigma_inv(-t,s));


    %find x locations where D large
    D_cut = .5;
    lambda_cut = 0.5;

    %parameter values
    D_large = 1;
    D_small = D_large*1e-2;

    lambda_large = 0.25;
    lambda_small = lambda_large*1e-2;

    %construct boundary points, interior
    xm_int = 1:mn*xn;

    x_bd_0 = 1:mn;
    x_bd_l = (xn-1)*mn+1:xn*mn;

    x_bd = union(x_bd_0,x_bd_l);

    m_bd_0 = 1:mn:mn*xn;
    m_bd_n = mn:mn:xn*mn;

    m_bd = union(m_bd_0,m_bd_n);

    bd = union(m_bd,x_bd);

    x_int = 1:mn*xn;
    x_int(x_bd) = [];

    m_int = 1:total;
    m_int(m_bd) = [];

    xm_int(bd) = [];

    %extra bound for sensors
    m_bd_1_int = 1:mn-2:(mn-2)*(xn-1)+1;
    m_bd_nm1_int = (mn-2):mn-2:(mn-2)*xn;

    %now define the velocity (in the m-direction) given the above variables
    V = zeros(mn-1,1);
    for i = 1:mn-1
        V(i) = g((m(i)+m(i+1))/2);
    end

    %construct vectors for east and west velocity.
    Ve = V(2:end);
    Vw = V(1:end-1);
    Vw_m1 = V(end);
    Ve_m0 = V(1);

    %measure of implicitness (1 = backward euler, 0 = forward euler, 1/2 =
    %crank-nicholson)
    theta = 0.5;


    Vec = Ve'*dt/dm;
    Vec = repmat(Vec,1,xn);
    Vec = Vec(:);

    Vwc = Vw'*dt/dm;
    Vwc = repmat(Vwc,1,xn);
    Vwc = Vwc(:);

    Vwc_m1 = Vw_m1'*dt/dm;
    Vwc_m1 = repmat(Vwc_m1,1,xn);
    Vwc_m1 = Vwc_m1(:);

    Vec_m0 = Ve_m0'*dt/dm;
    Vec_m0 = repmat(Vec_m0,1,xn);
    Vec_m0 = Vec_m0(:);


    % 
    Vec = @(t) f(s(t))*Vec;
    Vwc = @(t) f(s(t))*Vwc;
    Vwc_m1 = @(t) f(s(t))*Vwc_m1;
    Vec_m0 = @(t) f(s(t))*Vec_m0;



    %initial condition
    u0 = 10/3;
    % IC = u0*(X<=6).*(X>=2).*(M>=.2).*(M<=.4);
    IC = IC_1_d_m(M).*(X<=5);%.*exp(-M);
    IC1d = sum(IC)/mn;
    IC = IC(:);



    %sigma for flux limiters
    sigma = @(r) (r+abs(r))./(1+abs(r));


    %determine indices where lambda, D large or small. 

    D_m_large = x_int(M(x_int)>=D_cut);
    D_m_small = x_int(M(x_int)<D_cut);

    D_bd_0_large = x_bd_0(M(x_bd_0)>=D_cut);
    D_bd_0_small = x_bd_0(M(x_bd_0)<D_cut);

    D_bd_l_large = x_bd_l(M(x_bd_l)>=D_cut);
    D_bd_l_small = x_bd_l(M(x_bd_l)<D_cut);


    ind_total = 1:total;

    lambda_large_ind = M<=.5;
    lambda_small_ind = M>.5;

    lambda_large_ind = lambda_large_ind(:);
    lambda_small_ind = lambda_small_ind(:);


    %Define sparse matrix as a function for computation
    %need to define for when velocity is both positive and negative.

    A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[vw.*(-1+sw/2); ...
        (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],total,total);

    A_pos_m1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
        (-vw.*sw/2)],total,total);

    A_pos_m0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*(1-1*se/2); ve.*se/2],total,total);

    Dmat_large = D_large*dt/dx^2*sparse([D_m_large D_m_large D_m_large],[D_m_large-mn D_m_large+mn D_m_large],...
        [ones(1,2*length(D_m_large)) -2*ones(1,length(D_m_large))],total,total);

    Dmat_large_bd = D_large*dt/dx^2*sparse([D_bd_0_large D_bd_0_large D_bd_l_large D_bd_l_large],[D_bd_0_large D_bd_0_large+mn D_bd_l_large D_bd_l_large-mn],...
            2*[-ones(1,length(D_bd_0_large)) ones(1,length(D_bd_0_large)) -ones(1,length(D_bd_l_large))...
            ones(1,length(D_bd_l_large))],total,total);

    Dmat_small = D_small*dt/dx^2*sparse([D_m_small D_m_small D_m_small],...
        [D_m_small-mn D_m_small+mn D_m_small],[ones(1,2*length(D_m_small))...
        -2*ones(1,length(D_m_small))],total,total);

    Dmat_small_bd = D_small*dt/dx^2*sparse([D_bd_0_small D_bd_0_small D_bd_l_small D_bd_l_small],...
        [D_bd_0_small D_bd_0_small+mn D_bd_l_small D_bd_l_small-mn],...
        2*[-ones(1,length(D_bd_0_small)) ones(1,length(D_bd_0_small))...
        -ones(1,length(D_bd_l_small)) ones(1,length(D_bd_l_small))],total,total);

    Dmat = Dmat_large + Dmat_large_bd + Dmat_small + Dmat_small_bd; %incorporate both interior and boundary for simplicity

    A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-vw.*sw/2); ...
        (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],total,total);

    A_neg_m1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[(-vw.*sw/2); ...
        (vw.*sw/2-vw)],total,total);

    A_neg_m0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*se/2; (ve-ve.*se/2)],total,total);



    integ_ind = [];
    for i = 1:xn
        integ_ind = [integ_ind i*ones(1,mn)]; 
    end
    
    add_matrix = sparse(integ_ind,1:mn*xn,1); %matrix with 1's along rows for addition

    integ_matrix = add_matrix;

    %initialize
    u = zeros(total,tn);
    u(:,1) = IC;

    tic

    for i = 2:tn


        %Change sensors, A matrix based on velocity (whose sign depends on
        %f(s(t))
        if f(s(t(i))) >= 0 
            %sensors
            [r_e,r_w,r_w_m1,r_e_m1] = positive_sensor(u(:,i-1),m_int,m_bd_n,m_bd_0,m_bd_1_int,m_bd_nm1_int);
            %construct matrices
            A_int = A_pos(sigma(r_e),sigma(r_w),Vec(t(i)),Vwc(t(i)),m_int,1);
            A_m1 = A_pos_m1(sigma(r_w_m1),Vwc_m1(t(i)),m_bd_n,1);
            A_m0 = A_pos_m0(sigma(r_e_m1),Vec_m0(t(i)),m_bd_0,1);
        else
            %sensors
            [r_e,r_w,r_w_m1,r_e_m1] = negative_sensor(u(:,i-1),m_int,m_bd_n,m_bd_0,m_bd_1_int,m_bd_nm1_int);
            %construct matrices
            A_int = A_neg(sigma(r_e),sigma(r_w),Vec(t(i)),Vwc(t(i)),m_int,1);
            A_m1 = A_neg_m1(sigma(r_w_m1),Vwc_m1(t(i)),m_bd_n,1);
            A_m0 = A_neg_m0(sigma(r_e_m1),Vec_m0(t(i)),m_bd_0,1);
        end

        %integrate over m (just riemann for now)
    %     w = dm*accumarray(integ_mat,u(:,i-1));
        w = dm*integ_matrix*u(:,i-1);
        w = repmat(w',mn,1);
        w = w(:);

        %main computation
        u(:,i) = (speye(total) + theta*(A_int + A_m1 + A_m0))\((speye(total)...
            - (1-theta)*(A_int + A_m1 + A_m0) + Dmat)*u(:,i-1) + dt*(lambda_large*u(:,i-1).*lambda_large_ind +...
            lambda_small*u(:,i-1).*lambda_small_ind).*(1-w)); 

    end

    toc

    %y for visualization
    y = zeros(mn,xn,tn);
    z = zeros(xn,tn);

    for i = 1:tn
        y(:,:,i) = reshape(u(:,i),mn,xn);
        z(:,i) = dm*sum(y(:,:,i));
    end


end
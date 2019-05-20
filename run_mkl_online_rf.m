
function [w_new,f_hat,theta]=run_mkl_online_rf(y,D,x,w,params)
N=params.N; T=params.T; ker_list=params.ker_list;   beta=params.beta; delta=params.delta;
eta=params.eta;
lambda=params.lambda;



n_ker=length(params.sigma);
theta=params.theta_ini;
%w=cell(size(ker_list));

    
for i=1:n_ker


    w0=w{i};
    %ker=ker_list{i};
    %kx{i}=kernelmatrix(ker,D,x,1,0,2);
    %kx0=kx{i};
    kx0=[sin(D{i}*x); cos(D{i}*x)];
    fx(i)=kx0'*w0;
    er_temp(i,1)=(y-fx(i))^2;


   % theta(i)=theta(i)*beta^er_temp+0.0001;


    l_g=-2*(y-fx(i))*kx0;


    w_temp=(1-eta*lambda/10)*w0-eta*l_g/10;
    w0=w_temp;
    %kx{i}=[kx{i};kernelmatrix(ker,x,x,1,0,2); ];
   % fx_new(i)=kx0'*w0;
    w{i}=w0;
end
%theta=theta/sum(theta);
f_hat=fx*theta/sum(theta);
%ratio=er_temp/min(er_temp+0.001);
   mx=max(er_temp);
    mn=min(er_temp);
    mm=(mx+mn)/2;
    mr=(mx-mn)/2;
    
    ratio=(er_temp-mm+1e-12)/(mr+1e-12);
theta=theta.*beta.^ratio+1e-12;
%theta=theta/sum(theta);
    if max(theta)>1e10
    theta=theta/sum(theta);
    end

%D_new=[D x];
w_new=w;    
end



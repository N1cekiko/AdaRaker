
function [w_new,f_hat,theta]=run_mkl_online_rf(y,D,x,w,params)
N=params.N; T=params.T; ker_list=params.ker_list;   beta=params.beta; delta=params.delta;
eta=params.eta;
lambda=params.lambda;


n_ker=length(params.sigma);
theta=params.theta_ini;
%w=cell(size(ker_list));

    
for i=1:n_ker


    w0=w{i};

    kx0=[sin(D{i}*x); cos(D{i}*x)];
    fx(i)=(kx0'*w0);
    er_temp(i,1)=max(0,sign(-y*fx(i)));

    l_g=grad_log(y,fx(i))*kx0;


    w_temp=(1-eta*lambda)*w0-eta*l_g;
    w0=w_temp;
    w{i}=w0;
end

theta=theta/sum(theta);
f_hat=(fx*theta);
theta=theta.*(beta.^er_temp)+0.0001;

w_new=w;    
end



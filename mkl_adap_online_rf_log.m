

function [f_hat, theta_new,act_w,er_hat,erm]=mkl_adap_online_rf(y,X,params)
S=params.S;
N=params.N; T=params.T; ker_list=params.ker_list;   beta=params.beta; delta=params.delta;
lambda=params.lambda; sigma=params.sigma; L=params.L;
%lambda=params.lambda;
%D_new=X(:,1);
n_act=1; N_act=floor(log(T)/log(S))+1; n_ker=length(sigma); 
theta=params.theta_ini;
w=cell(size(ker_list));
W_int=zeros(2*L,1);
D=cell(n_ker,1); w_new=cell(N_act,1);F_hat=zeros(N_act,T);theta_new=cell(N_act,1);
eta=zeros(N_act,1);act_w=zeros(N_act,1);
for i=1:n_ker
 D{i}=sigma(i)*randn(L,N);
% w_new{i}=W_int;
end


C=ceil(L/N)-1;
resd=L-C*N;
for i=1:n_ker
    for c=1:C+1
        G=randn(N,N);
        [Q,R]=qr(G);
        %D_temp(c*N-N+1:c*N,:)=sigma(i)*Q;
        if c<C+1
            D_temp(c*N-N+1:c*N,:)=sigma(i)*Q;
        else
            D_temp(c*N-N+1:c*N-N+resd,:)=sigma(i)*Q(1:resd,:);     
        end
    end
    v1=chi2rnd(ones(L,1));
 D{i}=diag(sqrt(v1))*D_temp;
end

for k=1:N_act
    theta_new{k}=theta;
    w_new{k}=cell(n_ker,1);
    temp=cell(n_ker,1);
    for i=1:n_ker
        temp{i}=W_int;
    end
    w_new{k}=temp;
    act_w(k)=1;
end
for t=	1:T     
    yt=y(t);
    x=X(:,t);
    for k=1:n_act
        eta(k)=min(1/S,1/sqrt(S^k));
        if  mod(t-1,S^(k-1))==0 
            params.eta=min(1/S,1/sqrt(S^k));
            act_w(k)=params.eta;

            w=w_new{k};
                
            params.theta_ini=theta_new{k};
          
            
        else
            w=w_new{k};
            params.theta_ini=theta_new{k};
            params.eta=min(1/S,1/sqrt(S^k));
        end
        [w_temp,f_temp,theta_temp]=run_mkl_online_rf_log(yt,D,x,w,params);
         w_new{k}=w_temp;F_hat(k,t)=f_temp;theta_new{k}=theta_temp;
         er(k,1)=max(0,sign(-yt*f_temp));
      %  er(k,1)=abs(f_temp-yt);
    end
   

    f_hat(t)=sign(act_w'*F_hat(:,t)/sum(act_w));
    er_hat(t)=max(0,sign(-yt*f_hat(t)));
    if t==1
        er_hat(t)=1e0;
    end
    mx=max(er);
    mn=min(er);
    mm=(mx+mn)/2;
    mr=(mx-mn)/2;
    
    ratio=(er-mm+1e-10)/(mr+1e-10);
    beta=params.eta;
    act_w(1:n_act)=act_w(1:n_act).*beta.^ratio+1e-10;
%    if max(act_w)>1e12
        act_w=act_w/sum(act_w);
 %   end
 %   act_w=act_w/sum(act_w);
%     f_hat(t)=act_w'*F_hat(:,t)/sum(act_w);
%     er_hat(t)=(f_hat(t)-yt)^2;
%     if t==1
%         er_hat(t)=1;
%     end
    erm(t)=mean(er_hat);
      n_act=floor(log(t)/log(S))+1;
end
% figure(1)
% plot(erm)
% 
% figure(2)
% plot(gm)
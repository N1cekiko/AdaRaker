clear all;
close all;



load TomData


T=length(y);
N=size(X,1);
 

    all_beta=0.5;
    % all_lambda=st*T;
    all_lambda=0.01;
    all_cplx=50;
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Strongly adaptive OMKL with random feature approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
params5=struct;
params5.N=N; params5.T=T;  params5.ker_list={'rbf'}; n_ker=size(params5.ker_list,2);
params5.beta=all_beta; params5.delta=0.1;    params5.w_ini=0*ones(1,n_ker); 
params5.lambda=all_lambda;   params5.sigma=[ 0.1 1 10]; params5.L=all_cplx;  params5.theta_ini=ones(length(params5.sigma),1)/length(params5.sigma);
params5.S=2;


    tic;
    [f_rf,theta_rf,w_hat,er_rf,erm_rf]=mkl_adap_online_rf(y,X,params5);        
    t_adaprf=toc;
    E_adaraker=erm_rf(end);
    
    





    fprintf('AdaRaker=%f\n\n', E_adaraker);

    

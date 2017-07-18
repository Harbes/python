% ========================================================================
%
%      Program to estimate a linear equation by GMM
%
% ========================================================================
function gmm_opt(  )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );

    t = 1000;    % Choose the sample size      

    xt = rand(t,1);
    ut = randn(t,1);
    yt = 10 + 5*xt + ut;


    lmax = 2;    %  Length of the Newey-West Covariance matrix estimation 
    crit  = 0.00001;                      %* Convergence criterion         
    maxit = 30;                           %* Maximum number of iterations  

    % OLS estimates
    b0 = [ones(length(xt),1) xt]\yt;        
    disp(['OLS parameter estimates  = ', num2str(b0') ]);
    
    % First stage estimation 
    ops = optimset('LargeScale','off','Display','off');
    flag = 0;
    [ bgmm,qmin] = fminunc(@(b) gmmcrit(b,yt,xt,flag,lmax),b0,ops);
  
    % Second stage estimation
    flag = 1;
    [ bgmm,qmin] = fminunc(@(b) gmmcrit(b,yt,xt,flag,lmax),bgmm,ops);
    
    [obj,w] = gmmcrit(bgmm,yt,xt,flag,lmax);
    dg      = numgrad(@meaneqn,bgmm,yt,xt);    
    v       = dg'*inv(w)*dg;
    vinv    = inv(v)/t;

    disp(' ');
    disp(['The value of the objective function = ', num2str(obj) ]);
    disp(['J-test                              = ', num2str(t*obj) ]);
    disp(['GMM estimates                       = ', num2str(bgmm') ]);
    disp(['GMM standard errors                 = ', num2str(sqrt(diag(vinv))') ]);
    disp(['GMM t-statistics                    = ', num2str(bgmm'./sqrt(diag(vinv))') ]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Procedure which defines the moment equations  
%-------------------------------------------------------------------------
function mt = meqn(b,yt,xt)

    t  = length(yt);
    ut = yt - b(1) - b(2)*xt;                 
    zt = [ones(t,1),xt];   
        
    [~,cols] = size( zt );      
    k = 1;        
    mt = zeros(length(zt),2);
    for j = 1:cols
        
        mt(:,k) = zt(:,j).*ut;
        k=k+1;
    end  
                              

end
%-------------------------------------------------------------------------
% Procedure which defines the mean of the moment conditions   
%-------------------------------------------------------------------------
function ret = meaneqn(b,yt,xt)

        ret = mean(meqn(b,yt,xt))';

end
%-------------------------------------------------------------------------
% GMM objective function which also computes the optimal W
%-------------------------------------------------------------------------
function [ret,w] = gmmcrit(b,yt,xt,flag,lmax)
        
        d = meqn(b,yt,xt);
        t = length(d);
        g = mean(d)';

        if flag
           
            w  = d'*d;

            tau=1;

            while tau <= lmax
                wtau = d((tau+1):size(d),:)'*d(1:(size(d)-tau),:);
                w    = w + (1.0-tau/(lmax+1))*(wtau + wtau');
                tau  = tau + 1;
            end
            w = w./t;         
        else
             w = eye(size(d,2));
        end
        ret = g'*inv(w)*g;

end



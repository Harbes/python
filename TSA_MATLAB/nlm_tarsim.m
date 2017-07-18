%=========================================================================
%
%   Program to demonstrate the sampling properties of TAR models
%
%=========================================================================
function nlm_tarsim( )

    clear all;
    clc;

    % Initialise the random number generator  
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) ); 	

    % Choose model to simulate
    % 1 = SETAR, 2 = STAR, 3 = LSTAR, 4 = ESTAR
    type = 3; 
    
    % Choose test type
    % p = 1: test based on regressing yt on {constant ylag, ylag*ylag}       
    % p = 2: test based on regressing yt on {constant ylag, ylag*ylag, ylag*ylag^2}   
    % p = 3: test based on regressing yt on {constant ylag, ylag*ylag, ylag*ylag^2, ylag*ylag^3}  
    % p = 4: test based on regressing yt on {constant ylag, ylag*ylag,              ylag*ylag^3}  
    p = 3;

    % Parameters
    c      = 0.0;                                                     
    phi1   = -0.5;                                           
    beta1  = [0.0 0.5 1.0];      
    gam    = [0.5 1.0 1.5];                      
    sig    = 5;                        
    nobs   = 100;                                                     
    ndraws = 10000;                                       

    % Initialise arrays
    lm       = zeros(ndraws,1);
    lm_power = zeros(size(beta1,1),1);

    for j = 1:length(beta1)

        for k=1:ndraws
        
            %  Simulate a lstar model
            y           = tar(type,nobs,phi1,beta1(j),gam(1),c,sig);                
            [lm(k),dof] = tar_test(y,p);                                     

        end
                       
        % Size
        if j==1
            
            cve = quantile(lm,.95);
            cva = chi2inv(0.95,dof);
            pow = 100*mean(lm>cve);   
            nom = 100*mean(lm>cva);

            disp(' ')
            disp(['Size of test: beta1             = ', num2str(beta1(j)) ]);
            disp(['5%  critical value (empirical)  = ', num2str(cve) ]); 
            disp(['5%  critical value (asymptotic) = ', num2str(cva) ]);  
            disp(['Nominal size                    = ', num2str(nom) ]);           
            lm_power(j) = pow;

        else
            
            disp(' ')
            disp(['Power of test: beta1        = ', num2str(beta1(j)) ]);
            disp(['Unadjusted                  = ', num2str(100*mean(lm>chi2inv(0.95,dof)))]);
            disp(['Size adjusted               = ', num2str(100*mean(lm>cve)) ]);
            lm_power(j) = 100*mean(lm>cve)';

        end

    end

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Simulate a Threshold Autoregressive models
%-------------------------------------------------------------------------
function y = tar(gtype,nobs,phi1,beta1,gam,c,sig)

   
    y = zeros(nobs+100,1) ;
    u = randn(nobs+100,1);

    for  t = 2:nobs+100

        % SETAR
        if gtype == 1
            g = y(t-1) > c; 
        end                                                                 
        % STAR
        if gtype == 2
            g = normcdf( gam*(y(t-1) - c) ); 
        end                                                                   
        % LSTAR
        if gtype == 3
            g = 1/( 1 + exp(-gam*(y(t-1) - c)) ); 
        end                                                                  
        % ESTAR
        if gtype == 4
            g = 1 - exp(-gam*(y(t-1) - c)^2); 
        end                                                                  
        y(t) = phi1*y(t-1) + beta1*y(t-1)*g + sig*u(t);
    end
    y = trimr(y,100,0);
end


%-------------------------------------------------------------------------
%  Compute the LM statistic to test a threshold autoregressive model 
%  assuming one lag in the auxiliary model
%-------------------------------------------------------------------------
function [ lm,dof ] = tar_test(yvar,p)

    % First stage regression      
    y = trimr(yvar,1,0);
    x = [ones(size(y,1),1) , trimr(yvar,0,1)];                     
    k = size(x,2);
    u = y - x*(x\y);                                        

    % Second stage regression      
    if p == 1
        x = [x  x(:,2).*(x(:,2))]; 
    elseif p == 2
        x = [x  x(:,2).*(x(:,2)).^1  x(:,2).*(x(:,2).^2)]; 
    elseif p == 3
        x = [x , x(:,2).*(x(:,2).^1) , x(:,2).*x(:,2).^2 , x(:,2).*x(:,2).^3]; 
    elseif p == 4
        x = [x , x(:,2).*(x(:,2).^1) , x(:,2).*(x(:,2).^3)]; 
    end     
    e = u - x*(x\u);                                        

    % Compute LM statistic        
    r2  = 1 - sum(e.^2)/sum( (u-mean(u)).^2 );
    lm  = length(y)*r2;
    dof = size(x,2) - k;  
    
end





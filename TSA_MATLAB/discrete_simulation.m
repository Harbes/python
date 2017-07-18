%=========================================================================
%
%   Program to demonstrate alternative qualitative response models
%
%=========================================================================
function discrete_simulation( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );
    ops = optimset('LargeScale','off','Display','off'); 

    % Simulate data set 
    theta0 = [1 ; 1 ];
    t      = 100;
    ndraws = 5000;

    theta_full   = zeros(ndraws,2);
    theta_probit = zeros(ndraws,2);
    theta_tobit  = zeros(ndraws,2);
    theta_trunc  = zeros(ndraws,2);
    theta_rest   = zeros(ndraws,2);

    for j = 1:ndraws;

        u = randn(t,1);
        x = [ ones(t,1) randn(t,1) ];
        y = x*theta0 + u;

          
        % Probit data
        y_1 = ones(t,1);
        ind = y < 0.0;
        y_1(ind)  = 0.0;                    

        % Tobit data
        ind  = y > 0.0;
        y_2  = y.*ind;          

        % Truncated data
        y_3 = y(ind);
        x_3   = x(ind,:);

        % Full data set results  
        theta_full(j,:)= x\y;           

        % Probit data set results  
        theta_probit(j,:) = fminunc( @(b) lprobit(b,y_1,x),theta0,ops);
        
        % Tobit data set results  
        theta_tobit(j,:) = fminunc( @(b) ltobit(b,y_2,x),theta0,ops);
    
        % Truncated data set results   
        theta_trunc(j,:) = fminunc( @(b) ltrunc(b,y_3,x_3),theta0,ops);

        % Restricted data set results  
        theta_rest(j,:) = x_3\y_3;
    end

    disp(['True parameter values = ',num2str(theta0(1)),'  ',num2str(theta0(2))]);

    disp(' ');

    disp(['Mean (full)           = ',num2str(mean(theta_full))]);
    disp(['Mean (Probit)         = ',num2str(mean(theta_probit))]);
    disp(['Mean (Tobit)          = ',num2str(mean(theta_tobit))]);
    disp(['Mean (Trunctated      = ',num2str(mean(theta_trunc))]);
    disp(['Mean (Restricted)     = ',num2str(mean(theta_rest))]);

    disp(' ');
    
    
    disp(['Bias (full)           = ',num2str(mean(theta_full)-theta0')]);
    disp(['Bias (Probit)         = ',num2str(mean(theta_probit)-theta0')]);
    disp(['Bias (Tobit)          = ',num2str(mean(theta_tobit)-theta0')]);
    disp(['Bias (Truncated      = ',num2str(mean(theta_trunc)-theta0')]);
    disp(['Bias (Restricted)     = ',num2str(mean(theta_rest)-theta0')]);

    disp(' ');
    
    tmp = bsxfun(@minus,theta_full,theta0');
    disp(['RMSE (full)           = ',num2str(mean(tmp).^2)]);
    tmp = bsxfun(@minus,theta_probit,theta0');
    disp(['RMSE (Probit)         = ',num2str(mean(tmp).^2)]);
    tmp = bsxfun(@minus,theta_tobit,theta0');
    disp(['RMSE (Tobit)          = ',num2str(mean(tmp).^2)]);
    tmp = bsxfun(@minus,theta_trunc,theta0');
    disp(['RMSE (Truncated)      = ',num2str(mean(tmp).^2)]);
    tmp = bsxfun(@minus,theta_rest,theta0');
    disp(['RMSE (Restricted)     = ',num2str(mean(tmp).^2)]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Probit negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = lprobit(b,y,x)

        f  = normcdf(x*b);
        lf = -mean( y.*log(f) + (1 - y).*log(1 - f) );

end
%-------------------------------------------------------------------------
%  Tobit negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = ltobit(b,y,x)

        m = x*b;
        d = y ~= 0;
        
        lf = -mean( d.*log(normpdf(y - m)) + (1 - d).*log(1 - normcdf(m)) );

end
%-------------------------------------------------------------------------
%  Truncated negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = ltrunc(b,y,x)
     
     m  = x*b;
     lf = -mean( log(normpdf(y - m)) - log(1 - normcdf(0.0 - m)) );
     
end


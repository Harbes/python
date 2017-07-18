%=========================================================================
%
%    Simulating the power of the Wald test using exponential regression
%    model
%
%=========================================================================
function test_power( )

    clear all;
    clc;
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) )

    beta0       = 1;            % Intercept parameter                                             
    t           = 5;            % Sample size                                 
    ndraws      = 10000;
    wd          = zeros(ndraws,1);
    beta1_range = [-4, -3, -2, -1, 0, 1, 2, 3, 4];
    power       = zeros(1,length(beta1_range));

    for k = 1:length(beta1_range);
    
        beta1 = beta1_range(k);
    
        j  = 0; 
        x  = randn(t,1);      % Explanatory variable (fixed in repeated samplea)    

        for i = 1:ndraws          % Main do loop to generate Monte Carlo results   
            mue     = exp(beta0 + beta1*x);  %mean
            u       = rand(t,1);
            y       = -mue.*log(1 - u);
   
            theta0  = [beta0, beta1];
            options = optimset('LargeScale', 'off', 'Display', 'off');
            theta   = fminunc(@(theta) lnl(theta,y,x),theta0,options);
            H       = analytic_hessian(theta,y,x);
   
            % One restriction
            R       = [0 1];
            Q       = 0;
            wd(i,1) = t*(R*theta' - Q)'*inv( R*inv(-H)*R' )*(R*theta' - Q);

            if wd(i)>4.288
            
                j = j+1;
            end
    
        end
  
    power(k) = j/ndraws;

    end
    disp(['beta1 =   ' num2str(beta1_range)]);
    disp(['size for =' num2str(power)]);

end



%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
% Define the unconstrained log of the likelihood at each observations
%-------------------------------------------------------------------------

function ln1 = lnlt(theta,y,x)

     beta0 = theta(1);
     beta1 = theta(2);
     mue   = exp(beta0 + beta1*x);
     ln1   = -log(mue) - y./mue;
     
end
%-------------------------------------------------------------------------
% Define the unconstrained log of the likelihood 
%-------------------------------------------------------------------------
function ln = lnl(theta,y,x)

    lnlt1 = lnlt(theta,y,x);
    ln = -mean(lnlt1);

end


%-------------------------------------------------------------------------
% Define the analytic Hessian
%-------------------------------------------------------------------------

function H = analytic_hessian(theta,y,x)

    beta0 = theta(1);
    beta1 = theta(2);
    mue  = exp(beta0 + beta1*x);

    H      = zeros(2,2);
    H(1,1) = -sum(y./mue);
    H(1,2) = -sum((x.*y)./mue );
    H(2,1) = H(1,2);
    H(2,2) = -sum( ((x.^2).*y)./mue );
    
end
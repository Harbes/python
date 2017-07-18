% ========================================================================
%
%     Estimating a model of heteroskedasticity 
%
% ========================================================================
function hetero_estimate( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

    % Simulate the model   
    t       = 500;
    [y,x,w] = simulatedata( t );
                       
    % Estimate the unconstrained model  
    theta = [1; 2; 0.1; 0.1];
    [theta1,l1,~,~,~,H1] = fminunc(@(b) neglog1(b,y,x,w),theta);

    disp(['Log-likelihood (unconstrained) = ', num2str(-l1) ]);
    disp( 'Unconstrained parameter estimates ' );
    disp( [theta1 (1/t)*sqrt( diag(inv(H1) ) ) ] );
    
    % Estimate the constrained model      
    theta = [1; 2; 0.1 ];
    [theta0,l0,~,~,~,H0] = fminunc(@(b) neglog0(b,y,x,w),theta);

    disp(['Log-likelihood (constrained) = ', num2str(-l0) ]);
    disp( 'Constrained parameter estimates ' );
    disp( [theta0 (1/t)*sqrt( diag(inv(H0) ) ) ] );
 
end
%
% -------------------------- Functions ---------------------------------
%
%-----------------------------------------------------------------------
%   Simulate the data  
%-----------------------------------------------------------------------
function [ y,x,w ] = simulatedata( t )

    beta0 = 1;
    beta1 = 2;
    gam0  = 0.1;
    gam1  = 0.1;

    x = randn(t,1);                                               
    w = (0.0:0.1:0.1*(t-1))';                                                

    u = sqrt(exp(gam0 + gam1*w)).*randn(t,1);                        
    y = beta0 + beta1*x + u;                                

end
%-----------------------------------------------------------------------
% Negative unconstrained log-likelihood  
%-----------------------------------------------------------------------
function lf = neglog1(theta,y,x,w)

    lf = -mean( lnlt1(theta,y,x,w) );

end
%-----------------------------------------------------------------------
% Unconstrained log-likelihood function at each observation
%-----------------------------------------------------------------------
function lnl = lnlt1(b,y,x,w)

    mu  = b(1) + b(2)*x;
    sig = sqrt( exp(b(3) + b(4)*w) );
    lnl = -(1/2)*log(2*pi*sig.^2) - (y - mu).^2 ./(2*sig.^2);       


end
%-----------------------------------------------------------------------
% Negative constrained log-likelihood  
%-----------------------------------------------------------------------
function lf = neglog0(theta,y,x,w)

    lf = -mean( lnlt0(theta,y,x,w) );

end
%-----------------------------------------------------------------------
% Constrained log-likelihood function at each observation
%-----------------------------------------------------------------------
function lnl = lnlt0(b,y,x,w)

    mu   = b(1) + b(2)*x;
    sig  = sqrt( exp(b(3) + 0*w) );
    lnl  = -(1/2)*log(2*pi*sig.^2) - (y - mu).^2 ./(2*sig.^2);        

end


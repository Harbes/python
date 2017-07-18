% ========================================================================
%
%    Testing a model of heteroskedasticity 
%
% ========================================================================

function hetero_test( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

    % Simulate the model   
    t       = 500;
    [y,x,w] = simulatedata( t );
                       
    % Estimate the unconstrained model  
    theta = [1; 2; 0.1; 0.1];
    [theta1,l1,~,~,~,H1] = fminunc(@(b) neglog1(b,y,x,w),theta);

    l1 = -l1;
    disp(['Log-likelihood (unconstrained) = ', num2str(l1) ]);
    disp( 'Unconstrained parameter estimates ' );
    disp( [theta1 (1/t)*sqrt( diag(inv(H1) ) ) ] );
    
    % Estimate the constrained model      
    theta = [1; 2; 0.1 ];
    [theta0,l0,~,~,~,H0] = fminunc(@(b) neglog0(b,y,x,w),theta);

    l0 = -l0;
    disp(['Log-likelihood (constrained) = ', num2str(l0) ]);
    disp( 'Constrained parameter estimates ' );
    disp( [theta0 (1/t)*sqrt( diag(inv(H0) ) ) ] );

    % LR test   
    lr = -2*t*(l0 - l1);
    disp(['LR statistic            = ',num2str(lr) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lr,1)) ]);

    % Wald test   
    r = [0 , 0 , 0 , 1];
    q = 0;
    wd = t*(r*theta1 - q)'*inv(r*inv(H1)*r')*(r*theta1 - q);
    disp(['Wald statistic          = ',num2str(wd) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',wd,1)) ]);


    % LM test      
    theta = [theta0 ; 0];
    gmat  = numgrad(@lnlt1,theta,y,x,w);  
    g     = mean(gmat)';
    j     = gmat'*gmat/t;
    lm    = t*g'*inv(j)*g; 
    disp('Gradient evaluated at contrained estimates');
    disp( g );
	disp('Outer product of gradients matrix');
    disp( j );     
    disp(['LM statistic            = ',num2str(lm) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lm,1)) ]);                

    % LM test (regression form)   
    x = [ones(t,1) , x];                                 
    % Stage 1 regression
    b = x\y; 
    u = y - x*b;    
    w = [ones(t,1) , w];                               
    v = u.^2;
    % Stage 2 regression
    b  = w\v; 
    e  = v - w*b;
    r2 = 1 - sum(e.^2)/sum( (v-mean(v)).^2 );
    lm = t*r2;
    disp(['LM statistic (regression) = ',num2str(lm) ]); 
    disp(['p-value                   = ',num2str(1-cdf('chi2',lm,1)) ]);

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

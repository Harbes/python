%=========================================================================
%
%     Program to perform tests on a nonlinear consumption function
%     
%     U.S. data on real consumption and real disposable income(2005 $)
%     1960:Q1 to 2009:Q4 200 observations
%
%=========================================================================

function nls_contest( ) 		 

    clear all; 
    clc;
 
    % Load data
    load USdata.mat

    yt = inc;
    ct = cons;
    t  = length( yt );

    % Estimate the constrained model 
    b0 = [-228;  0.9] ;     
    [theta0,f0] = fminunc(@(b) neglog0(b,ct,yt),b0);  

    disp('Restricted Parameter Estimates');
    theta0 = [ theta0' 1.000];
    disp( theta0');
    
    % Estimate the unconstrained model
    b0 = [theta0 ];
    [theta1,f1] = fminsearch(@(b) neglog1(b,ct,yt),b0); 
    
    f0 = -f0;
    f1 = -f1;
    
    % Compute relevant matrices
    g    = numgrad(@lnlt1,theta0',ct,yt); 
    G    = mean( g );
    J    = g'*g/t;
    invH = inv(numhess(@neglog1,theta1',ct,yt));
    
    disp('Unestricted Parameter Estimates');
    disp( theta1 );
    disp('Covariance Matrix of the Parameters');
    disp( invH/t );

    % Perform likelihood ratio test     
    lr = -2*t*(f0 - f1);               
    disp('Restricted Likelihood Function');
    disp( t*f0 );
    disp('Unrestricted Likelihood Function');
    disp( t*f1 );
    disp(' ');
    disp('LR test and p-value');
    disp( [lr 1-cdf('chi2',lr,1)] );
    
    % Perform Wald test           
    R = [0 0 1];
    Q = 1;
    W = t*(R*theta1' - Q)'*inv(R*invH*R')*(R*theta1' - Q);
    disp(' ');
    disp('Wald statistic and p-value');
    disp( [ W  1-cdf('chi2',W,1) ] );

    % Perform Lagrange multiplier test (OPG)    
    disp('');
    disp( 'Outer product of gradients matrix')
    disp( J );   
    LM = t*(G*inv(J)*G');                         
    pv = 1-cdf('chi2',LM,1);                           
    disp(' ');
    disp('LM statistic and p-value');
    disp( [ LM  1-cdf('chi2',LM,1) ] );

    % Perform 2-step LM test    
    x = [ ones(t,1),inc ];
    y = cons;
    % Estimate constrained model
    b = x\y;         
    e = y - x*b;             
    b = [b ; 1.0];

    % Evaluate derivatives at constrained estimates
    z1 = ones(t,1);                           
    z2 = inc.^b(3);
    z3 = b(2)*log(inc).*inc.^b(3);
    z  = [z1 , z2 , z3];
    
    % Second step regression
    v  = e - z*(z\e);                                                   
    r2 = 1 - e'*e\v'*v;                                                  

    % LM statistic
    lm = t*r2;                 
    disp('2- step LM statistic and p-value');
    disp( [ lm  1-cdf('chi2',lm,1) ] );

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Log-likelihood function of constrained model
%-------------------------------------------------------------------------
function lf = neglog0(b,ct,yt)
     
     t  = length( ct );
     e  = ct - b(1) - b(2)*yt;                                 
     s2 = e'*e/t;                               % Concentrate the variance      
     f  = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e.^2/s2 ;
     lf = -mean( f );

end
%-------------------------------------------------------------------------
% Log-likelihood function of unconstrained model
%-------------------------------------------------------------------------
function lf = neglog1(b,ct,yt)
     
     t  = length(ct);
     e  = ct - b(1) - b(2)*yt.^b(3);                                
     s2 = e'*e/t;                              % Concentrate the variance   
     f  = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e.^2/s2;
     lf = -mean( f );
     
end
%-------------------------------------------------------------------------
% Log-likelihood function of constrained model at each observation
%-------------------------------------------------------------------------
function lf = lnlt1(b,ct,yt)
     
     t  = length(ct);
     e  = ct - b(1) - b(2)*yt.^b(3);                                
     s2 = e'*e/t;                              % Concentrate the variance   
     lf = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e.^2/s2;
     
end

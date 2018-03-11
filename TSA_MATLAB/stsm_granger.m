%=========================================================================
%
%   Program to perform Granger causality tests on a VAR
%
%=========================================================================
function stsm_granger( )

    clear all
    clc

    % Read the data: quarterly US data from Jan-1959 to Dec-1998
    load sims_data.mat
 
    % Define variables
    r    = ytdata(:,1);        
    lex  = log( ytdata(:,2) );
    lcp  = log( ytdata(:,3) );
    lm   = log( ytdata(:,4) );
    lp   = log( ytdata(:,5) );
    lo   = log( ytdata(:,6) );
    sdum = ytdata(:,7:17);
    
    t = length(ytdata);
    
    % Construct variables for use in VAR
    % interest rate and the annual percentage growth rates in money, price and output
    yvar = [ r   lm   lp   lo ];    
    tmp  = 100*(trimr(yvar(:,2:4),12,0) - trimr(yvar(:,2:4),0,12));
    y    = [ trimr(yvar(:,1),12,0) tmp]; 

    % Granger causality tests (based on the Wald test)
    t   = length(y);
    
    % Obtain estimates of the VAR(2) by OLS to each equation  
    bols   = [ones(t-2,1) trimr(y,1,1) trimr(y,0,2)]\trimr(y,2,0);      
    theta0 = bols(:);                                                   

    % Estimate model: will converge in one step
    ops = optimset('LargeScale','off','Display','iter');
    [theta1,fval,~,~,~,H] = fminunc(@(b) neglog( b,y ),theta0,ops);                                        

    cov1 = (1/t)*inv(H);
%%    
    % Wald test of no causality from money to interest rates

    r      = zeros(2,36);
    r(1,3) = 1.0;
    r(2,7) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof    = size(r,1);

    disp(' ')
    disp(['Wald test (money to interest rates)  = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from prices to interest rates
    r      = zeros(2,36);
    r(1,4) = 1.0;
    r(2,8) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
   
    disp(' ')
    disp(['Wald test (prices to interest rates) = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from output to interest rates
    r      = zeros(2,36);
    r(1,5) = 1.0;
    r(2,9) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
   
    disp(' ')
    disp(['Wald test (output to interest rates) = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from all to interest rates  
    r      = zeros(2,36);
    r(1,3) = 1.0;
    r(2,7) = 1.0;
    r(3,4) = 1.0;
    r(4,8) = 1.0;
    r(5,5) = 1.0;
    r(6,9) = 1.0;
    q      = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof    = size(r,1);
    disp(' ')
    disp(['Wald test (all to interest rates)    = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from interest rates to money
    r      = zeros(2,36);
    r(1,11) = 1.0;
    r(2,15) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (interest rates to money)  = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from prices to money
    r      = zeros(2,36);
    r(1,13) = 1.0;
    r(2,17) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (prices to money)          = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')
 
    % Wald test of no causality from output to money    
    r      = zeros(2,36);
    r(1,14) = 1.0;
    r(2,18) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (output to money)          = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    
    % Wald test of no causality from all to money    
    r      = zeros(2,36);
    r(1,11) = 1.0;
    r(2,15) = 1.0;
    r(3,13) = 1.0;
    r(4,17) = 1.0;
    r(5,14) = 1.0;  
    r(6,18) = 1.0;
    q      = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof    = size(r,1);
    disp(' ')
    disp(['Wald test (all to money)             = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from interest rates to prices  
    r      = zeros(2,36);
    r(1,20) = 1.0;
    r(2,24) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (interest rates to prices) = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')
    
    % Wald test of no causality from money to prices  
    r      = zeros(2,36);
    r(1,21) = 1.0;
    r(2,25) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (money to prices)          = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from output to prices  
    r      = zeros(2,36);
    r(1,23) = 1.0;
    r(2,27) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (output to prices)         = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from all to prices   
    r      = zeros(2,36);
    r(1,20) = 1.0;
    r(2,24) = 1.0;
    r(3,21) = 1.0;
    r(4,25) = 1.0;
    r(5,23) = 1.0;  
    r(6,27) = 1.0;
    q      = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof    = size(r,1);
    disp(' ')
    disp(['Wald test (all to prices)            = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from interest rates to output  
    r      = zeros(2,36);
    r(1,29) = 1.0;
    r(2,33) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (interest rates to output) = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % Wald test of no causality from money to output
    r       = zeros(2,36);
    r(1,30) = 1.0;
    r(2,34) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (money to output)          = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')
 
    % Wald test of no causality from prices to output   **/
    r      = zeros(2,36);
    r(1,31) = 1.0;
    r(2,35) = 1.0;
    q      = [ 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof = size(r,2);
    
    disp(' ')
    disp(['Wald test (prices to output)         = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')
 
    % Wald test of no causality from all to output
    r       = zeros(2,36);
    r(1,29) = 1.0;
    r(2,33) = 1.0;
    r(3,30) = 1.0;
    r(4,34) = 1.0;
    r(5,31) = 1.0;  
    r(6,35) = 1.0;
    q      = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
    wd     = (r*theta1 - q)'*inv(r*cov1*r')*(r*theta1 - q);
    dof    = size(r,1);
    disp(' ')
    disp(['Wald test (all to output)            = ',num2str(wd) ]);
    disp(['Number of degrees of freedom         = ',num2str(dof) ]);
    disp(['p-value                              = ',num2str(1-chi2cdf(wd,dof)) ]);

end

%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for an unrestricted VARMA model
%--------------------------------------------------------------------------
function f = neglog( b,y )

    [ t,n ] = size( y );
    v      = zeros( t,4 );
    lf      = zeros( t-2,1 );
    
    % First loop over MA part  
    for i = 3:t
        v(i,1) = y(i,1) -  b(1) -  b(2)*y(i-1,1) -  b(3)*y(i-1,2) -  ...
            b(4)*y(i-1,3) -  b(5)*y(i-1,4) -  b(6)*y(i-2,1) -        ...
            b(7)*y(i-2,2)  - b(8)*y(i-2,3)  - b(9)*y(i-2,4);

        v(i,2) = y(i,2) - b(10) - b(11)*y(i-1,1) - b(12)*y(i-1,2) - ...
            b(13)*y(i-1,3) - b(14)*y(i-1,4) - b(15)*y(i-2,1) -      ...
            b(16)*y(i-2,2) - b(17)*y(i-2,3) - b(18)*y(i-2,4);

        v(i,3) = y(i,3) - b(19) - b(20)*y(i-1,1) - b(21)*y(i-1,2) - ...
            b(22)*y(i-1,3) - b(23)*y(i-1,4) - b(24)*y(i-2,1) -      ...
            b(25)*y(i-2,2) - b(26)*y(i-2,3) - b(27)*y(i-2,4);

        v(i,4) = y(i,4) - b(28) - b(29)*y(i-1,1) - b(30)*y(i-1,2) - ...
            b(31)*y(i-1,3) - b(32)*y(i-1,4) - b(33)*y(i-2,1) -      ...
            b(34)*y(i-2,2) - b(35)*y(i-2,3) - b(36)*y(i-2,4);
    end
    v  = trimr( v,2,0 );    
    vc = v'*v/length(v);

    for i = 1:t-2;

        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(vc)) ...
                - 0.5*v(i,:)*inv(vc)*v(i,:)';    
    end
    f = -mean( lf ); 
end

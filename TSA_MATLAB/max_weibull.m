%=========================================================================
%
%     Program to estimate the parameters of a Weibull distribution using
%     the Newton-Raphson and BHHH algorithms
%
%=========================================================================

function max_weibull( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) )

    t = 20;

    % Generate Weibull random numbers
    alpha = 1.0;
    beta  = 2.0;
 
    z = rand(t,1); 
    x = -log(1 - z)/alpha;  % Generate exponential random numbers  
    y = x.^(1/beta);        % Generate Weibull random numbers 

    % or use Matlab function 
    y = wblrnd((1/alpha)^(1/beta),beta,[t,1]);  % Generate Weibull random numbers

    % or load data
    y = [0.293, 0.589, 1.374, 0.954, 0.608, 1.199, 1.464, ...
        0.383, 1.743, 0.022, 0.719, 0.949, 1.888, 0.754, ...
        0.873, 0.515, 1.049, 1.506, 1.090, 1.644]';
    
    disp('      t        y');
    disp([(1:1:t)',y]);

    % Starting values 
    theta0 = [0.5;1.5];             
    alpha  = theta0(1);
    beta   = theta0(2);

    g = gradvec(theta0,y);
    h = hessmat(theta0,y);
    j = opgmat(theta0,y);

    disp( 'Estimates at iteration 0 ' );
    disp( theta0 );
    disp(['Value of log-likelihood at iteration 0    = ' num2str(-lnl(theta0,y))] );
    disp('Gradient at iteration 0' ) 
    disp( g );
    disp('Hessian at iteration 0');
    disp( h );
    disp('OPG at iteration 0');
    disp(j);

    % Newton Raphson update
    theta1 = theta0 - inv(h)*g;    

    disp( 'Estimates at iteration 1 (Newton Raphson) ' );
    disp( theta1 );
    disp(['Value of log-likelihood at iteration 1    = ' num2str(-lnl(theta1,y))] );
    
    % BHHH
    theta1 = theta0 + inv(j)*g;           

    disp( 'Estimates at iteration 1 (BHHH) ' ); 
    disp( theta1 );
    disp(['Value of log-likelihood at iteration 1    = ' num2str(-lnl(theta1,y))] );

    % Call fminunc to optimize function
    theta = fminunc(@(theta) lnl(theta,y),theta1);
    
    disp( 'MLE of theta ' );
    disp( theta );

    g = gradvec(theta,y);
    h = hessmat(theta,y);
    j = opgmat(theta,y);

    disp( 'Covariance matrix (Hessian) ' ); 
    disp( (1/t)*inv(h) );

    disp( 'Covariance matrix (OPG) ' );
    disp( (1/t)*inv(j) );


    % Call fminunc to optimize the concentrated likelihood
    beta = fminunc(@(b) lnlc(b,y),theta1(2));
    
    alpha = 1/mean(y.^beta);
    disp( 'Results for concentrated likelihood' );
    disp(['MLE of alpha   = ' num2str(alpha) ] );
    disp(['MLE of beta    = ' num2str(beta) ] );
    
    % Call fminunc to optimize the transformed likelihood
    [theta,~,~,~,~,h] = fminunc(@(theta) lnlt(theta,y),theta1);

    disp( 'Results for transformed likelihood' );
    disp(['MLE of lambda   = ' num2str(theta(1)) ] );
    disp(['MLE of beta     = ' num2str(theta(2)) ] );
    
    % Std error by delta method
    d = [ -(1/theta(2))*theta(1)^(-1/theta(2)-1), ...
        (log(theta(1))/theta(2)^2)*theta(1)^(-1/theta(2)) ];   
    disp( 'Standard error of lambda by delta method' ) ;
    disp( -d*inv(h)*d' );
end

%------------------------- Functions ----------------------------- %

% Log-likelihood 
function lf =  lnl(theta,y)
 
     a  = theta(1);
     b  = theta(2);
     f  = log(a) + log(b) + (b-1)*log(y) - a*y.^b;
     lf = -mean( f );

end

% Concentrated log-likelihood 
function lf =  lnlc(b,y)
 
     a  = 1/mean(y.^b);    
     f  = log(a) + log(b) + (b-1)*log(y) - a*y.^b;
     lf = -mean( f );

end

% Transformed log-likelihood function
function lf = lnlt(theta,y)

     l  = theta(1);
     b  = theta(2);
     f  = log(b) - log(l ) + (b-1)*log(y/l) - (y/l).^b;
     lf = -mean(f);
     
end

% Return gradient vector 
function g = gradvec(theta,y)

    alpha = theta(1);
    beta  = theta(2);

    g1 = 1/alpha - mean(y.^beta);
    g2 = 1/beta + mean(log(y)) - alpha*mean(log(y).*y.^beta );
    g  = [g1; g2];             

end

% Return Hessian matrix 
function h = hessmat(theta,y)

    alpha = theta(1);
    beta  = theta(2);

    h11 = - 1/alpha^2 ;
    h12 = - mean(log(y).*y.^beta);
    h21 = h12;
    h22 = -1/beta^2 - alpha*mean( (log(y).^2).*(y.^beta) );
    h   = [h11, h12; h21, h22];
end

    % Return OPG matrix 
function j  = opgmat(theta,y)

    alpha = theta(1);
    beta  = theta(2);

    g22 = 1/beta + log(y) - alpha*log(y).*y.^beta;
    j11 = mean( ( 1/alpha - y.^beta).^2 );
    j12 = mean( ( 1/alpha - y.^beta) .* g22 );
    j21 = j12;
    j22 = mean( g22.^2 );
    j = [j11, j12; j21, j22];
    
end


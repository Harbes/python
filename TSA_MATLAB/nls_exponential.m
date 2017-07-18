%==========================================================================
%
%    Estimate the exponential model by maximum likelihood.
%     
%==========================================================================

function nls_exponential() 		 

    clear all;
    clc;
    
    t = 50;
    %flag = 1;       % 1 = simulate data
    flag = 0;        % 0 = use GAUSS data to reproduce text
   
    if flag 
        [ y,x ] = simulatedata( t );
    else      
        % Aternatively load the GAUSS data
        load nls_expdata.dat
        y = nls_expdata(:,1);
        x = nls_expdata(:,2);     
    end

    % Estimate the model and compute Hessian se  
    start = [0.1,0.1];   
    theta = fminunc(@(theta) neglog(theta,y,x),start);      

    ht   = numhess(@neglog,theta',y,x);
    vcov = (1/t)*inv(ht);
    
	disp(['Parameter estimates = ',num2str(theta) ]);
	disp('Negative Hessian matrix');
    disp( ht );
    disp('Covariance matrix');
    disp( vcov );
end

%--------------------- Functions ----------------------%
% Simulate data
function [ y,x ] = simulatedata( t )

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );
    
    b0  = 1.0;
    b1  = 0.05;
    sig = 0.5;

    u = sig*randn(t,1);

    x = [1:1:t]';

    y = b0*exp( b1*x ) + u;
    
end

% Log-likelihood function
function lf = neglog(b,y,x)

    lf = - mean( lnlt(b,y,x) );
    
end
% Log-likelihood (concentrated) at each observation
           
function lf = lnlt(b,y,x)

     e  = y - b(1)*exp( b(2)*x );                              
     s2 = e'*e/length(e);                        
     lf = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e.^2/s2 ;

end

%--------------------------------------------------------------------------
%   Define the log of the likelihood (unconcentrated) at each observation
%--------------------------------------------------------------------------
       

function val = lnlt_unc(b,y,x,t)

     e  = y - b(1)*exp( b(2)*x );   % Residual error                          
     s2 = b(3);                     % Do not concentrate out the residual variance   
     val =  - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e.^2/s2 ;

end
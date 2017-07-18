%=========================================================================
%
%   Program to estimate a stochastic frontier model 
%
%=========================================================================

function nls_frontier(  )

    clear all;
    clc;
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) )

    % Set parameter values
    beta0 = 1;         
    beta1 = 0.5;        
    sig1  = 1;          
    sig2  = 1.5;          

    t = 1000;     
    x = randn(t,1); 
    
    % Compute density and cumulative distribution function
    h = 0.01;        
    u = seqa(-10,h,1501);
    f = (1/sig2)*exp( sig1^2/(2*sig2^2) + u/sig2 ).*normcdf( -(u + sig1^2/sig2)/sig1 );

    %********************************************************************
    %***
    %***     Generate graph
    %***
    %********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    plot(u,f,'-k')
    ylabel('$f(u)$');
    xlabel('$u$');
    set(gca,'YTick',[0.0 0.1 0.2 0.3]');
    
    box off    
   
    % Print the tex file to the relevant directory
    %laprint(1,'frontier','options','factory');
        
    % Monte Carlo simulation
    ndraws = 5000;
    theta  = zeros(ndraws,4);
    
    % cdf to be used as a look-up table in inverse cdf method 
    cdf    = h*cumsum(f);                                                 

    y = zeros(t,1);
    for i = 1:ndraws
        
        % Generate random numbers for y  
        for j = 1:t
            
            [~,ind] = min(abs(cdf - rand));
            y(j)       = beta0 + beta1*x(j) + u(ind);   
            
        end
         
        % Estimate the model
        theta0     = [beta0 ; beta1 ; sig1 ; sig2 ];      
        bhat       = fminsearch(@(b) neglog(b,y,x),theta0);
        theta(i,:) = abs(bhat);
        
    end


disp( '                            beta0      beta1      sig1       sig2');
disp( ['Population parameters = ', num2str(theta0') ]);
disp( ['Mean                  = ', num2str(mean(theta)) ]);
%print "Bias                  = " meanc(theta_mle)'-theta0';
%print "MSE                   = " meanc((theta_mle-theta0').^2)';



    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%   Log-likelihood function of the stochastic frontier model
%-------------------------------------------------------------------------
function lnf = neglog(b,y,x) 

    u    = y - b(1) - b(2)*x;
    sig1 = abs(b(3));
    sig2 = abs(b(4));
        
    % Normal - exponential likelihood
    lf = - log(sig2) + 0.5*sig1^2/sig2^2 + u/sig2 ...
           + log(normcdf( -(u + sig1^2/sig2 )/sig1 ) );

    % Normal - halfnormal likelihood
%        sigs = sqrt( sig2^2 + sig1^2 );
%        lam  = sig2/sig1;
%        lf   = 0.5*log(2/pi) - log(sigs) + log(cdfn( -s*e*lam/sigs )) ...
%               - 0.5*e.^2/sigs^2 ;       

    lnf = -mean( lf );

end



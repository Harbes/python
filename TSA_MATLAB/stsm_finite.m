%=========================================================================
%
%   Simulation demonstration of the finite sample properties of the 
%   AR(1) estimator
%
%=========================================================================

function stsm_finite( )

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234567) );

    % Generate the data
    t      = 50;
    sigv   = 1.0;
    phi1   = 0.9;
    ndraws = 10000;
    
    theta_mle = zeros( ndraws,1 );
    
    for i = 1:ndraws
        
        vt = sqrt(sigv)*randn( t+101,1 );
        yt = zeros( length( vt ),1 );
        
        % Simulate the AR(1) model
        for j = 2:length( vt )
            
            yt(j) = phi1*yt(j-1) + vt(j);
        end
        
        % Get rid of first 100 observations
        yt = trimr( yt,100,0);  
        
        % Conditional mle
        theta_mle(i,:) = trimr( yt,0,1 )\trimr( yt,1,0 );
        
    end
    
    % Compute statistics of sampling distribution
    mse_mle  = mean( (theta_mle - phi1).^2 );     
    rmse_mle = sqrt( mse_mle );
        
    disp( ['Population parameter    = ', num2str( phi1 ) ] )
    disp( ['Mean (cond. mle)        = ', num2str( mean(theta_mle) ) ] )
    disp( ['Bias (cond. mle)        = ', num2str( mean( theta_mle )- phi1 ) ] )
    disp( ['Bias (Shenton-Johnson)  = ', num2str( -2*phi1/t ) ] )
    disp( ['RMSE (cond. mle)        = ', num2str( rmse_mle ) ] )
    
    
end
        
  
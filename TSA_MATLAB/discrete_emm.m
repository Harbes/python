%=========================================================================
%
%   Sampling distribution of EMM estimator for the INARMA(1,1) model
%   derived using Monte Carlo methods.
%
%   THIS TAKES AN AGE TO RUN ... reduce niter, ndraws &/or n to get an idea
%
%=========================================================================
function discrete_emm( )

    clear all
    clc
    
    RandStream.setGlobalStream( RandStream('mt19937ar','seed',123) );
    
    t      = 50;                                                     
    rho    = 0.3;                                                                  
    beta   = 0.7;
    lam    = 3.5;        
    theta0 = [rho ; lam ; beta ];    
    
    % Control length of simulation run
    h      = 50;     
    n      = t*h;             

    % Number of replications to generate finite sample properties
    ndraws = 5000;
    
    % Maximum number of searches
    niter  = 200;                                                      

    % Loop over auxiliary models and lag length
    for m = 1:2

        % Start lag at k=2 for identification 
        for k = 2:3 
            

            model = m;
            lag   = k;

            % Main DO LOOP to generate sampling distribution
            theta_emm  = zeros(ndraws,3);

            for j =1:ndraws 
                
                disp(['# draws  = ',num2str(j) ])


                % Generate the actual data for binomial thinning model 
                % Generate Poisson random variables
                u = poissrnd(lam,t+100,1);
                % Initialize y using the unconditional distribution 
                y = poissrnd(lam*(1+beta)/(1-rho),t+100,1);    
                z = y;                                                                                                                     
                for i = 2:t+100 
                    % AR part
                    if z(i-1) >0;                           

                        % Generate t_t-1 uniform random numbers    
                        v = rand(z(i-1),1);           
                        % Sum bernoulli random draws  
                        b_thin = sum( v < rho );      
                        % Generate realisation of y at time t   
                        z(i) = b_thin + u(i);           
                    else
                        % There are no survivors
                        z(i) = u(i);                    
                    end
                
                    % MA part
                    if u(i) >0  
                        % Generate uniform random numbers    
                        v = rand(u(i),1);           
                        % Sum bernoulli random draws  
                        b_thin = sum( v < beta );      
                        % Generate realisation of y at time t   
                        y(i) = b_thin + z(i-1);           
                    else
                        % There are no survivors
                        y(i) = z(i-1);                    
                    end
                end
            
                % Trim first 100 observations to overcome startup problems
                y = trimr(y,100,0);                    
            
                % Conditional least squares
                x    = [ ones(length(y)-1,1) trimr(y,0,1) ];
                bhat = x\trimr(y,1,0);
                ehat = trimr(y,1,0) - x*bhat;
                s2   = mean(ehat.^2);

                theta_0 = [ bhat(2) ; bhat(1) ; 0.1 ];            
    
                % Estimate the auxiliary model using actual data  
                if lag == 1;

                    x    = [ ones(length(y)-lag,1) trimr(y,0,1) ];
                    bhat = x\trimr(y,1,0);
                    ehat = trimr(y,1,0) - x*bhat;
                    s2   = mean(ehat.^2);

                elseif lag == 2;

                    x    = [ones(length(y)-lag,1)   trimr(y,1,1)   trimr(y,0,2) ];
                    bhat = x\trimr(y,2,0);
                    ehat = trimr(y,2,0) - x*bhat;
                    s2   = mean(ehat.^2);
            
                elseif lag == 3;

                    x    = [ones(length(y)-lag,1)  trimr(y,2,1)  trimr(y,1,2)  trimr(y,0,3)];
                    bhat = x\trimr(y,3,0);
                    ehat = trimr(y,3,0) - x*bhat;
                    s2   = meanc(ehat.^2);
            
                end

                % Choose auxiliary model type  
                if model == 1
                    g = bsxfun(@times,ehat,x);                
                elseif model == 2
                    g = [bsxfun(@times,ehat,x)  (ehat.^2 - s2) ]; 
                end

                % Compute the optimal weighting matrix      
                i = g'*g;
                p = 0.0;

                for l=1:p 

                    gam = g((l+1):length(g),:)'*g(1:(length(g)-l),:);		
                    i   = i + (1.0 - l/(p+1))*(gam + gam');				
                end
            
                i = i/length(g);
                iinv = inv(i);		

                [theta,~] = search(abs(theta_0),n,lag,model,iinv,t,niter,bhat);

                theta_emm(j,:) = theta;

            end


            % Sampling distribution of the emm estimator   
            mm    = mean(theta_emm);
            stdev = sqrt(bsxfun(@minus,theta_emm,mm).^2);
            rmse  = sqrt(bsxfun(@minus,theta_emm,theta0').^2);

            disp( ' ' )
            disp('EMM estimator results') 
            disp('------------------------------------------------')
            disp(['Number of replications              = ',num2str(ndraws) ]);
            disp(['Sample size                         = ',num2str(t) ]);
            disp(['Auxiliarly model                    = ',num2str(model) ]);
            disp(['Lags                                = ',num2str(lag) ]);
            disp(['Number of searches                  = ',num2str(niter) ]);
            disp(['Number of simulation paths          = ',num2str(h) ]);
            disp('True  Mean  Bias Std Dev RMSE')
            disp( [theta0 mm' mean(stdev)' mean(rmse)'])

        end

    end
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  The objective function to compute the EMM estimator
%-------------------------------------------------------------------------
function lf = q(b,n,lag,model,iinv,t,bhat)

        RandStream.setGlobalStream( RandStream('mt19937ar','seed',1) );

        rho_emm  = b(1);
        lam_emm  = b(2);
        beta_emm = b(3);

        % Generate simulated data    
        us = poissrnd(lam_emm,n,1);                  
        ys = poissrnd(lam_emm*(1+beta_emm)/(1-rho_emm),n,1);    

        zs = ys;                                                             
 
        for i = 2:n
            
            % AR part
            if zs(i-1) >0;                           

                % Generate t_t-1 uniform random numbers    
                vs = rand(zs(i-1),1);           
                % Sum bernoulli random draws  
                b_thin = sum( vs < rho_emm );      
                % Generate realisation of y at time t   
                zs(i) = b_thin + us(i);           
            else
                % There are no survivors
                zs(i) = us(i);                    
            end
            
            % MA part
            if us(i) >0  
                % Generate uniform random numbers    
                vs = rand(us(i),1);           
                % Sum bernoulli random draws  
                b_thin = sum( vs < beta_emm );      
                % Generate realisation of y at time t   
                ys(i) = b_thin + zs(i-1);           
             else
                    % There are no survivors
                    ys(i) = zs(i-1);                    
             end

        end

        if lag == 1

            xs    = [ones(length(ys)-lag,1)   trimr(ys,0,1) ];
            ehats = trimr(ys,1,0) - xs*bhat;

        elseif lag == 2

            xs    = [ ones(length(ys)-lag,1)   trimr(ys,1,1)  trimr(ys,0,2) ];
            ehats = trimr(ys,2,0) - xs*bhat;

        elseif lag == 3

            xs    = [ ones(length(ys)-lag,1)   trimr(ys,2,1)   trimr(ys,1,2)   trimr(ys,0,3)];
            ehats = trimr(ys,3,0) - xs*bhat;

        end
                    
        % Choose auxiliary model type      

        if model == 1
            gs = bsxfun(@times,ehats,xs);   
        elseif model == 2
            gs = [ bsxfun(@times,ehats,xs)   (ehats.^2 - s2) ]; 
        end
        
        f  = mean(gs)*iinv*mean(gs)';
        lf = -0.5*t*f;

end
%-------------------------------------------------------------------------
%  Random search function
%-------------------------------------------------------------------------
function [b0,f0] = search(b0,n,lag,model,iinv,t,niter,bhat)

    
    % Function evaluation at initial parameters
    f0 = q(b0,n,lag,model,iinv,t,bhat);                           

    b = zeros(length(b0),1);

    for i = 1:niter

        b(1) = rand;
        b(2) = rand*5;
        b(3) = rand;                              

        f = q(b,n,lag,model,iinv,t,bhat);

        if f < f0
            f0 = f; 
            b0 = b; 
        end
    end

end


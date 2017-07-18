%=========================================================================
%
%    Simulation example to generate the power function of the Wald test
%    of heteroskedasticity.
%                                                         
%=========================================================================

function hetero_power(  )
    
    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

    t = 50;  
    n = 10000;            
    x = rand(t,1);              

    % Population parameters                                    
    beta0  = 1.0;
    beta1  = 2.0;
    gam0   = 1.0;
    gam1   = 0.0:0.5:5;                      

    % Pre-allocate memory for estimates
    wd        = zeros(n,1);
    wd_power  = zeros( length(gam1),1);


    % Optimisation settings
    ops = optimset( 'LargeScale',  'off',    ...
                    'Display',     'off', ...
                    'MaxIter',     20000,    ...
                    'MaxFunEvals', 40000,    ...
                    'TolX',        5e-4,     ...
                    'TolFun',      5e-4          );
      

    % Restriction matrices
    R = [0  0  0  1 ];
    Q = 0;

    % Loop over elements of gamma
    for j = 1:length(gam1)
    
        start = [ beta0; beta1; gam0; gam1(j) ];
        % Simulation loop
        for k = 1:n

            % Generate data    
            u = sqrt( exp(gam0 + gam1(j)*x) ).*randn(t,1);   
            y = beta0 + beta1*x + u;                    

            % Estimate parameters
            [theta,~,~,~,~,H] = fminunc( @(p) neglog(p,y,x),start,ops );
        
            % Wald test
            wd(k) = t*(R*theta - Q)'*inv(R*inv(H)*R')*(R*theta - Q);

        end;
       
        if j==1; 
            
            c_5         = quantile(wd,.95);
            wd_h0       = wd; 
            wd_power(j) = 100*mean(wd > c_5 );

            disp('Size of test');
            disp('******************************');
            disp('5%  critical value (empirical)');
            disp(c_5);
            disp('Power ( 5% nominal sig. level)');
            disp( 100*mean( wd > chi2inv(0.95,1) ) );
            disp('Power ( 5% nominal sig. level)');
            disp( wd_power(j) );
            
        else
            
            wd_power(j) = 100*mean(wd > c_5 );
            disp('Power of test for gamma value');
            disp(gam1(j));
            disp('******************************');
            disp('Unadjusted power');
            disp( 100*mean( wd > chi2inv(0.95,1) ) );
            disp('Size adjusted power');
            disp(wd_power(j));

        end

    end
%**************************************************************************
%**
%**     Generate graphs
%**
%**************************************************************************
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    % Plot the sampling distribution of the Wald statistic under the null hypothesis
    hist(wd_h0,21)
    xlabel('Values of the Wald test statistic');
    ylabel('Empirical Distribution');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight;
    box off;

    %laprint(1,'waldsize','options','factory');

    % Plot the power function of the Wald statistic
    figure(2);

    plot(gam1,wd_power,'-k')
    xlabel('Values of the parameter $\gamma_1$');
    ylabel('Power (\%)');
    axis tight;
    box off;
    %laprint(2,'waldpower','options','factory');

end
%
%-------------------------Functions------------------------------------
%
%-----------------------------------------------------------------------
%      Negative log-likelihood function   
%-----------------------------------------------------------------------

function lf = neglog( p,y,x )

    lf = -mean( loglt(p,y,x) );
    
end
%-----------------------------------------------------------------------
%      Log-likelihood function at each observation
%-----------------------------------------------------------------------
function lft = loglt(p,y,x)

    mu  = p(1) + p(2)*x;
    sig = sqrt( exp(p(3) + p(4)*x) );
    lft = -(0.5)*log(2*pi*sig.^2) - ((y - mu).^2)./(2*sig.^2);        

end



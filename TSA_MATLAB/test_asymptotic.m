%==========================================================================
%
%   Program to generate asymptotic distribution of the Wald test
%   applied to the regression model.
%
%==========================================================================


function test_asymptotic(  )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

    % Set parameter values
    beta0  = 1.0;                                           
    beta1  = 0.0;                                                     
    beta2  = 0.0;                                                       
    beta3  = 0.0;                                                         
    sig2   = 0.1;                 
    sig    = sqrt(sig2);

    t      = 1000;                                                                     
    ndraws = 10000;

    x1 = rand(t,1);                                                                           
    x2 = randn(t,1);                                                                           
    x3 = randn(t,1).^2;                                                                        

    % Arrays to hold results
    wd1 = zeros(ndraws,1);                                                 
    wd2 = zeros(ndraws,1);                                                 
    wd3 = zeros(ndraws,1);                                               

    
    % Loop over number of replications
    theta0  = [beta0 ; beta1 ; beta2 ; beta3 ; sig2 ];
    options = optimset('LargeScale', 'off', 'Display', 'off');
   
    for i = 1:ndraws
       
        u = sig*randn(t,1);                                  
        y = beta0 + beta1*x1 + beta2*x2  + beta3*x3 + u;                               
        
        [p] =  fminsearch(@(p) neglog(p,y,x1,x2,x3),theta0,options); 
         
        H     = numhess(@neglog,p,y,x1,x2,x3);
        theta = p;
        cov   = inv(H);
        
        % One restriction
        R = [0  1  0  0  0]; 
        Q = 0;
        wd1(i) = t*(R*theta - Q)'*inv(R*cov*R')*(R*theta - Q);

        % Two restrictions
        R = [ 0  1  0  0  0 ;
              0  0  1  0  0 ]; 
        Q = [0 ; 0 ];
        wd2(i) = t*(R*theta - Q)'*inv(R*cov*R')*(R*theta - Q);
        
        % Three restrictions
        R = [ 0  1  0  0  0 ;
              0  0  1  0  0 ; 
              0  0  0  1  0 ]; 
        Q = [ 0 ; 0 ; 0 ];
        wd3(i) = t*(R*theta - Q)'*inv(R*cov*R')*(R*theta - Q);

    end
   
    %********************************************************************
    %***
    %***     Generate graph
    %***
    %********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
   
    [fcdf,x] = ecdf( wd3 ); 
    [f,bins] = ecdfhist(fcdf,x,21); 
    bar(bins,f,'hist');
    %[n,xout]=hist(rt,51);
    %bar(xout,n/t)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k');
    axis([0,15,-Inf,Inf])
    box off
    hold on
    ygrid = 0.0001:0.1:15;
    plot(ygrid,chi2pdf(ygrid,3),'-k','LineWidth',0.75)
    ylabel('$f(W)$');
    xlabel('$W$');
    %set(gca,'YTick',[] );
    hold off
    
    % Print the tex file to the relevant directory
    %laprint(1,'simwald','options','factory');

   
end

%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
% Log-likelihood function of unconstrained model
%-------------------------------------------------------------------------

function lf = neglog(theta,y,x1,x2,x3)

    m  = theta(1) + theta(2)*x1 + theta(3)*x2 + theta(4)*x3;    
    s2 = abs(theta(5));                                              
    lf = -mean( -0.5*log(2*pi) - 0.5*log(s2) - 0.5*(y - m).^2/s2  ) ;

end

% %------------------------------------------------------------------------- 
% %   Computes finite difference Hessian
% %------------------------------------------------------------------------- 
% 
% function H = numhess( f,x,varargin )
%     
%     k  = length( x );
%     f0 = feval( f, x, varargin{:} );
%  
%     % Compute the stepsize (h)
%     dx  = eps.^( 1/3 )*( abs(x) + eps );
%     xh  = x + dx;
%     dx  = xh - x;
%     ee  = diag( dx ); 
%  
%     % Compute forward and backward steps
%     fplus  = zeros( k,1 );
%     
%     for i=1:k
%         
%         fplus(i)  = feval( f, x+ee(:,i), varargin{:} );
%     end  
%     
%     H = zeros( k );
%     for j = 1:k
%        
%         for l = 1:j;
%             
%             H(j,l) = feval( f, x+ee(:,j)+ee(:,l), varargin{:} ); 
%         end
%     end
%     H = H + tril( H,-1 )';    
%     
%     fpp = bsxfun( @plus, fplus, fplus' ); 
%     H   = H - fpp + f0;
%     H   = H./( dx*dx' );    
%     
% end
% 

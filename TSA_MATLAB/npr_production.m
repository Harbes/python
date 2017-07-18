%=========================================================================
%
%    Gaussian nonparametric kernel regression estimator of a
%    production function using a product kernel. 
%
%    Example taken from Pagan and Ullah (1999) Non-Parametric
%    Econometrics.
%    The variables are yt (log of output), x1t and x2t (logs of inputs)
%
%=========================================================================

function npr_production( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

    % Simulate the nonlinear production function
    n = 2;                  % Number of explanatory variables
    t = 200;                % Number of observations

    xt = rand(t,n);        % Uniform random numbers  
    ut = 0.1*randn(t,1);
    yt = -0.2*log( exp(-5*xt(:,1)) + 2*exp(-5*xt(:,2)) ) + ut;

    disp('Standard deviation of x1t and x2t');
    disp(std(xt)); 

    % Compute the Nadaraya-Watson bivariate product kernel estimator
    xmesh   = (0.0:0.02:1)';
    [x1,x2] = meshgrid(xmesh);

    % True conditional mean
    mx_pop = -0.2*log( exp(-5*x1) + 2*exp(-5*x2) );    

    % Estimated conditional mean
    [mx,h] = kernbiv(yt,xt,x1(:),x2(:),0);  

    % Compute the derivatives of the kernel conditional mean at x1 = x2 = 0.5 
    
    % Compute correct bandwidths with s=1
    
    s=1;
    fac = -1/(4 + n + 2*s);
    h   = std(xt).*(t.^fac);

    [mx1,tmp] = kernbiv(yt,xt,0.5-h(1),0.5-h(2),1); 
    disp(mx1);
    [mx1,tmp] = kernbiv(yt,xt,0.5,0.5-h(2),1); 
    disp(mx1);
    [mx1,tmp] = kernbiv(yt,xt,0.5+h(1),0.5-h(2),1); 
    disp(mx1);

    [mx1,tmp] = kernbiv(yt,xt,0.5-h(1),0.5,1); 
    disp(mx1);
    [mx1,tmp] = kernbiv(yt,xt,0.5,0.5,1); 
    disp(mx1);
    [mx1,tmp] = kernbiv(yt,xt,0.5+h(1),0.5,1); 
    disp(mx1);

    [mx1,tmp] = kernbiv(yt,xt,0.5-h(1),0.5+h(2),1); 
    disp(mx1);
    [mx1,tmp] = kernbiv(yt,xt,0.5,0.5+h(2),1); 
    disp(mx1);
    [mx1,tmp] = kernbiv(yt,xt,0.5+h(1),0.5+h(2),1);
    disp(mx1);


    % x1 derivative
    [mxa,tmp] = kernbiv(yt,xt,0.5-h(1),0.5,1);
    [mxb,tmp] = kernbiv(yt,xt,0.5+h(1),0.5,1);
    d1 = (mxb-mxa)/(2*h(1));                          

    % x2 derivative
    [mxa,tmp] = kernbiv(yt,xt,0.5,0.5-h(2),1);
    [mxb,tmp] = kernbiv(yt,xt,0.5,0.5+h(2),1);
    d2 = (mxb-mxa)/(2*h(2));                          

    % Compute the derivatives of the population conditional mean
    xm1 = 0.5;      
    xm2 = 0.5;

    d1_pop = exp(-5*xm1)./(exp(-5*xm1) + 2*exp(-5*xm2));
    d2_pop = 2*exp(-5*xm2)./(exp(-5*xm1) + 2*exp(-5*xm2));


    fprintf(1,'\nPopulation derivative of x1 (evaluated at x1=0.5)    = %f\n',d1_pop); 
    fprintf(1,'Nonparametric derivative of x1 (evaluated at x1=0.5) = %f\n\n',d1); 
    
    fprintf(1,'Population derivative of x2 (evaluated at x2=0.5)    = %f\n', d2_pop); 
    fprintf(1,'Nonparametric derivative of x2 (evaluated at x2=0.5) = %f\n', d2); 
    

%**************************************************************************
%**
%**     Generate graphs
%**
%**************************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    %--------------------------------------------------------%
    
    mesh( x1,x2,reshape( mx,length(x1),length(x2) ), ...
      'EdgeColor','black');
    xlabel('$x_{1,t}$');
    ylabel('$x_{2,t}$');
    zlabel('$m(x_{1,t},x_{2,t})$');
    set(gca,'XTick',[0.0 0.2 0.4 0.6 0.8 1.0]);
    set(gca,'YTick',[0.0 0.2 0.4 0.6 0.8 1.0]);
    set(gca,'Ztick',[0.0 0.5 1.0])
    %axis tight
    grid 'off'
    box 'off'

    %laprint(1,'figprodn','options','factory');

end 
%
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%    Procedure to compute the Nadayara-Watson bivariate product kernel 
%    production function using a product kernel. 
%
%    Outputs are:
%
%     mx is a (n1 x n2) matrix of the conditional mean
%     h1 is the bandwidth for x1
%     h2 is the bandwidth for x2
%
%    Inputs are:
%
%     yt is a (t x 1) vector of the dependent variable
%     xt is a (t x 2) vector of the two independent variables
%     x1 is a (n1 x 1) vector of grid points associated with x1
%     x2 is a (n2 x 1) vector of grid points associated with x2
%
%-------------------------------------------------------------------------
function [ mx,h ] = kernbiv(yt,xt,x1,x2,s) 


    [t,n] = size(xt);
    
    x = [x1 x2];
    [tx,nx] = size(x);
    
    fac = -1/(4 + n + 2*s);
    h   = std(xt).*(t.^fac);
    ph  = prod( h );

    % Estimate density using product kernel
    ker  = zeros(n,1);
    pker = zeros(t,n);
    fx   = zeros(tx,1);
    fxy  = zeros(tx,1);
    for j = 1:tx
    
        for i = 1:t;
        
            for p = 1:n
                ker(p) = normpdf( (x(j,p) - xt(i,p))/h(p) );
            end 
            pker(i,1)  = prod( ker );
            pker(i,2) = prod( ker ).*yt(i);
        end
        fx(j)  = mean( pker(:,1) )/ph;
        fxy(j) = mean( pker(:,2) )/ph;
    end

    mx  = fxy./fx; 


end









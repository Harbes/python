%==========================================================================
%
%   	Program to estimate a symmetric BEKK MGARCH model of US yields
%       US daily zero coupon yields are expressed in percentages, 
%       starting 10/04/1988 and ending 12/28/2001
%
%==========================================================================
function mgarch_bekk( )

    clear all
    clc

    % Load data
    load yields_us

    % Choose variables
    r  = rdata(:,[ 1 2 ] );                                                  
    y  = 100*(trimr( r,1,0 ) - trimr( r,0,1 )); 
    t  = length( y );
 
    % Estimate the BEKK Aymmetric MGARCH(1,1) model             
    start = [
                1.0711572152073034 
                0.2847761383915719 
                0.5666734919888730 
                0.3849272367534938 
                0.0079982965205830 
                0.0492176759790856 
                0.2235088706616157 
                0.8810887055460768 
               -0.0064483299253031 
                0.0271203983383643 
                0.9730134169147020 
                0.0603486110799331 
                0.0587667747700864   ];
    ops = optimset( 'LargeScale','off','Display','iter' );
 
    [thetaa,lfa ] = fminunc( @(b) negloga(b,y),start,ops); 
    
    clear start
    lfa = -lfa;
    
    disp(['Likelihood function (asymmetric)   = ',num2str(lfa)]);
     
    % Estimate the BEKK Symmetric MGARCH(1,1) model  
    start = [
                1.0827864786474668 
                0.5341723472962869 
                0.6326817767441521 
                0.4018538161228361 
                0.0241620281075959
                0.2103118107333996 
                0.8982242699069298 
                0.0048863543542688 
                0.9609684149405032 
                0.0713182051952838 
                0.0474958290682067  ];
    ops = optimset( 'LargeScale','off','Display','iter' );
 
    [thetas,lfs ] = fminunc( @(b) neglogs(b,y),start,ops); 
    
    clear start
    lfs = -lfs;
    
    disp(['Likelihood function (symmetric)   = ',num2str(lfs)]);
    
    % Compute conditional variance at optimal parameters
    h = cov(y);     
    u = std(y);     
    u = u';
    
    c  = [ thetas(1)    0.0  ;     
           thetas(2)    thetas(3) ];

    a =  [ thetas(4)    thetas(5) ; 
           thetas(5)    thetas(6) ];

    d =  [ thetas(7)    thetas(8) ;
           thetas(8)    thetas(9) ];

    m = thetas([10 11]);                   
    
    cov_mat = zeros( t,4 );
    for i = 1:t
        
        u = y(i,:)'- m;                  % Update residuals    
        h = c*c' + a*(u*u')*a' + d*h*d'; % Update covariance matrix        
       
        cov_mat(i,:) = h(:)';
  
    end
      
    
    
    %*********************************************************************
    %**     Generate graph of conditional variances
    %*********************************************************************
    % Get daily dates from excel file
    %[num,txt] = xlsread('yields_us.xlsx','a2:g3308');
    %clear num
    %dvec     = datenum( txt(:,1), 'dd/mm/yyyy' );
    %vec      = dvec(2:end);
    vec = 1:t;
    
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    %--------------------------------------------------------%
    % Panel (a)
    subplot(2,2,1)
    plot( vec,cov_mat(:,1),'-k','LineWidth',0.75);
    %datetick('x','yyyy'); 
    title(' (a) Variance of 3-month rate' );
    %ylabel('$\widehat{h}_t$');
    xlabel('$t$');
    axis tight
    box off

      
    %--------------------------------------------------------%
    % Panel (b)
    subplot(2,2,2)
    plot( vec,cov_mat(:,4),'-k','LineWidth',0.75);
    %datetick('x','yyyy');
    title(' (b) Variance of 1-year rate' );
    %ylabel('$\widehat{h}_t$');
    xlabel('$t$');
    axis tight
    box off


    %--------------------------------------------------------%
    % Panel (c)
    subplot(2,2,3)
    plot( vec,cov_mat(:,2),'-k','LineWidth',0.75);
    %datetick('x','yyyy');
    title(' (c) Conditional covariance ' );
    %ylabel('$\widehat{h}_t$');
    xlabel('$t$');
    axis tight
    box off



    %--------------------------------------------------------%
    ccor = cov_mat(:,2)./( sqrt( cov_mat(:,1).*cov_mat(:,4) ) );
    % Panel (d)
    subplot(2,2,4)
    plot( vec,ccor,'-k','LineWidth',0.75);
    %datetick('x','yyyy');
    title(' (d) Conditional correlation ' );
    %ylabel('$\widehat{h}_t$');
    xlabel('$t$');
    axis tight
    box off
 
    %laprint(1,'bekk','options','factory');

    % Estimate the constant cov BEKK  MGARCH(1,1) model  
    start = [       1.215169555478360
                    0.781175687461714
                    0.400847056654352
                    0.401280406108247
                    0.213852259080924
                    0.899426993215549
                    0.963984125877018
                    0.073850907078739
                    0.040752440808561 ];
    
    ops = optimset( 'LargeScale','off','Display','iter' ); 
    [thetac,lfc ] = fminunc( @(b) neglogc(b,y),start,ops); 
    
    lfc = -lfc;
    
    % LR test of symmetry    
    lr  = -2*t*(lfs - lfa);

	disp(' ');
    disp(['LR test of symmetry   = ',num2str(lr)]);
	disp(['p-value               = ',num2str(1-chi2cdf(lr,4))]);
    
    % LR test of constant covariance
    lr  = -2*t*(lfc - lfa);

	disp(' ');
    disp(['LR test of const. cov.   = ',num2str(lr)]);
	disp(['p-value                  = ',num2str(1-chi2cdf(lr,4))]);

end
%
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-liklihood function for asymmetric BEKK uncontrained
%--------------------------------------------------------------------------
function lf = negloga( b,y )

    [ t,n ] = size( y );
    f      = zeros( t,1 );

    h = cov(y);     % Initialise  conditional covariance matrix
    u = std(y);     % Initialise the disturbance vector
    u = u';
    
    c  = [ b(1)    0.0  ;     
           b(2)    b(3) ];

    a =  [ b(4)    b(6) ; 
           b(5)    b(7) ];

    d =  [ b(8)    b(10) ;
           b(9)    b(11) ];

    for i = 1:t
        
        m    = b([12 13]);                   % Update conditional mean
        u    = y(i,:)'- m;                   % Update residuals                                
        f(i) = -0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*u'*inv(h)*u;        
        h    = c*c' + a*(u*u')*a' + d*h*d';  % Update conditional covariance matrix

    end
    
    lf = -mean( f );
    
end
%--------------------------------------------------------------------------
% Log-liklihood function for symmetric BEKK model
%--------------------------------------------------------------------------
function lf = neglogs( b,y )

    [ t,n ] = size( y );
    f      = zeros( t,1 );

    h = cov(y);     % Initialise  conditional covariance matrix
    u = std(y);     % Initialise the disturbance vector
    u = u';
    
    c  = [ b(1)    0.0  ;     
           b(2)    b(3) ];

    a =  [ b(4)    b(5) ; 
           b(5)    b(6) ];

    d =  [ b(7)    b(8) ;
           b(8)    b(9) ];

    for i = 1:t
        
        m    = b([10 11]);                   % Update conditional mean
        u    = y(i,:)'- m;                   % Update residuals                                
        f(i) = -0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*u'*inv(h)*u;       
        h    = c*c' + a*(u*u')*a' + d*h*d';  % Update conditional covariance matrix 


    end
    lf = -mean( f );
    
end
%--------------------------------------------------------------------------
% Log-liklihood function for constant covariance BEKK model
%--------------------------------------------------------------------------
function lf = neglogc( b,y )

    [ t,n ] = size( y );
    f      = zeros( t,1 );

    h = cov(y);     % Initialise  conditional covariance matrix
    u = std(y);     % Initialise the disturbance vector
    u = u';
    
    c  = [ b(1)    0.0  ;     
           b(2)    b(3) ];

    a =  [ b(4)    0.00 ; 
           0.00     b(5) ];

    d =  [ b(6)    0.00 ;
           0.00    b(7) ];

    for i = 1:t
        
        m    = b([8 9]);                   % Update conditional mean
        u    = y(i,:)'- m;                   % Update residuals                                
        f(i) = -0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*u'*inv(h)*u;       
        h    = c*c' + a*(u*u')*a' + d*h*d';  % Update conditional covariance matrix 


    end
    lf = -mean( f );
    
end



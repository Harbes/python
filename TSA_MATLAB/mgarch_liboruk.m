%==========================================================================
%
%   Program to estimate MGARCH model of UK libor rates based on an
%   SVAR model with GARCH conditional variances on the structural
%   disturbances and identification based on short-run restrictions.
%
%==========================================================================
function mgarch_liboruk( )
    
    clear all
    clc
    
    % Load data: 1 January 2004 to 28 May 2008
    % Missing observations in credit risk series fixed by interpoloation
    % 11/06/2004, 9/02/2005 & 26/12/2005
    
%     libor      = xlsread('libor_data','b3:e1152');
%     default    = xlsread('libor_data','f3:i1152');
%     repo       = xlsread('libor_data','j3:j1152');
%     swap       = xlsread('libor_data','k3:l1152');
%     liquidity  = xlsread('libor_data','m3:m1152');
%     volatility = xlsread('libor_data','n3:o1152');
    load libor_data.mat

    y = [  volatility      ...
           liquidity       ...
           swap            ...
           repo            ...
           default*10000   ... 
           libor ] ;
       
    % Dummy variable to allow for structural break in volatility 2/07/2007  
    dum = zeros( length(y),1 );
    dum(913:end) = 1;
    
    % Set SVAR parameters
    p = 2;      % Order of VAR lags
    q = 21;     % Order of VMA lags      

    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1);
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);
    s      = chol(omegav)';


   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);           


    bstart = [
                0.001855486206481
                0.068693371469487
                0.023880592273368
                0.001989809184491
                0.013648672676092
                0.000281803828058
               -0.121999669995451
               -0.002888991832387
                0.002650363820455
                0.217906436715095
                0.001195294344330
                0.008613191798914
                0.175400639477836
                1.103000246244063
                0.440791629588056
                0.495605700916872
                0.365097267442524
               -0.011572135693424
                0.650892298805934
               -1.176100342250688
               -0.342100219212525
               -1.059396487450620
               -0.545902063123713
               -0.729900334997059
               -0.402000678739595
                1.013901399027622
               -1.601300249119193
                1.018606365787122
                0.791195906556810
                0.818998254746848
                0.338395871264844
            ];

  options = optimset(   'Display',              'iter', ...
                        'MaxIter',              20000,   ...
                        'MaxFunEvals',          20000,   ...
                        'TolX',                 1e-8,   ...
                        'TolFun',               1e-4  );       
            
    [ b,fval ]  = fminunc( @(b) neglogl(b,v),bstart,options );     
  
    disp(' Log-likelihood function ');
    disp( fval );
    disp( b );
    
    % Matrix required to compute time-varying variance-covariances              
    s_store = getSmatrix( b,v );
    
    % Construct C(L) matrices (vector moving average)  
    c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ]);

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end
       tmp = ss';
       c  = [ c  tmp(:) ];

    end

    % Construct orthogonal impulse response functions
    tdate   = 1;          %Choose a point in time 
    s       = reshapeg(s_store(tdate,:),k,k)';
    impulse = vec(s');

    for i=2:q 

       tmp      = ( reshapeg( c(:,i),k,k )*s )';
       impulse  = [ impulse tmp(:) ];

    end 
	impulse = impulse';
    
       
    % Construct variance decompositions 
    tmp0   = reshapeg( cumsum( impulse.^2 ),q*k,k );
    tmp1   = repmat( sum(tmp0,2),1,k );
    decomp = reshapeg( 100*tmp0 ./ tmp1 , q , k^2 );
 
    format short
    disp('Variance decomposition (%) over time of libor in terms of factors');
    disp('*****************************************************************');
	disp(' ' );
    disp('Period         y1         y2        y3        y4         y5       y6 ');   

    disp( [ (1:1:size(decomp,1))' decomp(:,k*k-(k-1):k*k) ] );
    
    % Graph the instantaneous effect of shocks on the libor over time  
    rr            = size(v,1);
    decomp_time1  = zeros( rr,k );           % Lag 1 (instantaneous)      
    decomp_time5  = zeros( rr,k );           % Lag 5       
    decomp_time20 = zeros( rr,k );           % Lag 20

    for j = 1:rr
 
        s       = reshapeg(s_store(j,:),k,k)';
        impulse = vec(s');

        for i = 2:q

            tmp      = ( reshapeg( c(:,i),k,k )*s )';
            impulse  = [ impulse tmp(:) ];

        end 
        impulse = impulse';    
        tmp0    = reshapeg( cumsum( impulse.^2 ),q*k,k );
        tmp1    = repmat( sum(tmp0,2),1,k );
        decomp  = reshapeg( 100*tmp0 ./ tmp1 , q , k^2 );

        decomp_time1(j,:) = decomp(1,k*k-(k-1):k*k); 
        decomp_time5(j,:) = decomp(5,k*k-(k-1):k*k); 
        decomp_time20(j,:) = decomp(20,k*k-(k-1):k*k); 

    end    
    
    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    t = (2004:1/249:2008.61)';

    %--------------------------------------------------------%
    % Panel (a)
    subplot(3,3,1)
    plot(t,decomp_time1(:,1),'-k','LineWidth',1);
    title('Global Risk Factor');
    ylabel('One Period Lag');
    %ylim( [0 100] )
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );
    
    %--------------------------------------------------------%
    % Panel (b)
    subplot(3,3,2)
    plot(t,decomp_time1(:,3),'-k','LineWidth',1);
    title('Broad Liquidity Factor');
    %ylabel('Oil price');
    %ylim( [0 100] )
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );


    %--------------------------------------------------------%
    % Panel (c)
    subplot(3,3,3)
    plot(t,decomp_time1(:,6),'-k','LineWidth',1);
    title('Idiosyncratic Factor');
    %ylabel('Oil price');
    %ylim( [0 100] )
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );


    %--------------------------------------------------------%
    % Panel (d)
    subplot(3,3,4)
    plot(t,decomp_time5(:,1),'-k','LineWidth',1);
    %title('Money Shock');
    ylabel('Five Period Lag');
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );

    %--------------------------------------------------------%
    % Panel (e)
    subplot(3,3,5)
    plot(t,decomp_time5(:,3),'-k','LineWidth',1);
    %ylabel('Output');
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );


    %--------------------------------------------------------%
    % Panel (f)
    subplot(3,3,6)
    plot(t,decomp_time5(:,6),'-k','LineWidth',1);
    %ylabel('Output');
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );

    %--------------------------------------------------------%
    % Panel (g)
    subplot(3,3,7)
    plot(t,decomp_time20(:,1),'-k','LineWidth',1);
    ylabel('Twenty Period Lag');
    xlabel('t');
     box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );

    %--------------------------------------------------------%
    % Panel (h)
    subplot(3,3,8)
    plot(t,decomp_time20(:,3),'-k','LineWidth',1);
    %ylabel('Output');
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );

    %--------------------------------------------------------%
    % Panel (i)
    subplot(3,3,9)
    plot(t,decomp_time20(:,6),'-k','LineWidth',1);
    ylabel('Price');
    xlabel('t');
    box off
    axis tight
    set(gca,'YLimMode','manual')
    set(gca,'Ylim',[0 100]);
    set(gca,'YTick',[0 20 40 60 80 100] );
    set(gca,'XTick',[2004 2006 2008] );
    
    %--------------------------------------------------------%
    % Print the tex file to the relevant directory
    laprint(1,'decomp','options','factory');

end

%
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Wrapper function for log-liklihood
%--------------------------------------------------------------------------
function lf = neglogl( b,v )

    lf = - mean( loglt( b,v ) );
end

%--------------------------------------------------------------------------
% Log-likelihood function for SVAR
%--------------------------------------------------------------------------

function lf = loglt( b,v )


    % Unpack parameter vector
    delta = b(14:19).^2;
    alpha = normcdf( b(20:25) );
    beta  = normcdf( b(26:31) );
    
    [ t,n ] = size( v );
    lf      = zeros( t,1 );
    
    
    binv  = [   1        0      0     0       0     0  
                b(1)     1     b(8)   0       0     0  
                b(2)     0      1     0       0     0  
                b(3)     0      0     1       0     0  
                b(4)   b(6)   b(9)   b(11)    1     0  
                b(5)   b(7)   b(10)  b(12)  b(13)   1  ];

    u = ( inv(binv)*v' )';
            

    for i = 1:t;
        
        if i == 1;
            
            d = delta + alpha.*(std(u).^2)' + beta.*(std(u).^2)';
        else
            
            d = delta + alpha.*(u(i-1,:).^2)' + beta.*d;
        end
        
        omega = binv*diag(d)*binv';
        
        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(omega)) - 0.5*v(i,:)*inv(omega)*v(i,:)';
    end
    
    
end
%--------------------------------------------------------------------------
%   Return matrix to compute time-varying covariances
%--------------------------------------------------------------------------
function S = getSmatrix( b,v )


    % Unpack parameter vector
    delta = b(14:19).^2;
    alpha = normcdf( b(20:25) );
    beta  = normcdf( b(26:31) );
    
    [ t,n ] = size( v );
    S       = zeros( t,n^2 );
    
    
    binv  = [   1        0      0     0       0     0  
                b(1)     1     b(8)   0       0     0  
                b(2)     0      1     0       0     0  
                b(3)     0      0     1       0     0  
                b(4)   b(6)   b(9)   b(11)    1     0  
                b(5)   b(7)   b(10)  b(12)  b(13)   1  ];

    u = (inv(binv)*v')';
            
    for i = 1:t;
        
        if i == 1;
            
            d = delta + alpha.*(std(u).^2)' + beta.*(std(u).^2)';
        else
            
            d = delta + alpha.*(u(i-1,:).^2)' + beta.*d;
        end
       
        tmp    = binv*diag( sqrt(d) ); 
        
        S(i,:) = vec( tmp )';
              
    end
    
end

%--------------------------------------------------------------------------
% Reshape a matrix to agree with GAUSS reshape command
%--------------------------------------------------------------------------
function X = reshapeg(Y,r,c)

         tmp = reshape(Y',c,r);
         X   = tmp';

end



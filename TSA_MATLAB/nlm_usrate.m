%=========================================================================
%
%   Estimate an LSTAR model of US unemployment rate based on the 
%   simple specification to be found in Skalin and Terasvirta 
%
%=========================================================================
function nlm_usrate( )

    clear all
    clc
    
    % Load data 
    load usunemp.mat
    
    % Set up dates for plotting purposes
    startDate = datenum('01-31-1948');
    endDate   = datenum('03-31-2010');
    xData     = linspace(startDate,endDate,747);

%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    plot(xData,data,'-k','LineWidth',0.75);
    datetick('x','yyyy')
    ylabel('Unemployment Rate');
    axis tight
    box off

    %laprint(1,'unrate','options','factory');
      
    % Generate LSTAR variables    
    dy      = trimr(data,1,0) - trimr(data,0,1);    % dy(t)
    y_lag   = trimr(data,0,1);                      % y(t-1)
    dy12    = trimr(data,12,0) - trimr(data,0,12);  % y(t)-y(t-12)
    dy12lag = trimr(dy12,0,1);
    
    % Adjust variables to have same length as dy12lag
    dy    = trimr(dy,12,0);
    y_lag = trimr(y_lag,12,0);

    % Create data matrix 
    t = length(dy);     
    
    % Estimate the linear model
    xlin = [ ones( t,1 ) y_lag ];
    bols = xlin\dy;
    lr_linear = -bols(1)/bols(2);
    
    % Create elements for the nonlinear model
%     [ r,c ]       = size( xlin );
%     xnln          = (repmat( x(:,end),1,c ).^3).*xlin;
%     [ b,bint,u1 ] = regress( u0,[ xlin xnln ] );
%     RSS1          = sum( u1.^2);
% %     
% %     % Test for nonlinearity
%     test  = t*( (RSS0-RSS1)/RSS1);
%     disp( [test 1-chi2cdf( test,c )] );
%     
    % Optimization
    x       = [ dy y_lag dy12lag ];
    pstart  = 0.1*ones(6,1);
    
    options = optimset( 'Display',             'iter', ...
                        'MaxIter',              2000,   ...
                        'MaxFunEvals',          4000 );   
                    
    [phat,~,~,~,~,hess] = fminunc( @(p) logl(p,x),pstart,options );        
    
    vc = (1/t)*inv(hess);
    se = sqrt(diag(vc));
    
    disp(' Estimates    se       t-stat' );
    disp( [phat se phat./se] );
    
    lr_low = -phat(1)/phat(2);
    lr_high = -(phat(1)+phat(3))/(phat(2)+phat(4));

    disp('Long-run mean unemployment in linear state');
    disp( lr_linear );
    disp('Long-run mean unemployment in low state');
    disp( lr_low );
    disp('Long-run mean unemployment in high state');
    disp( lr_high ); 
     
end

%=========================================================================
%
%   Wrapper function
%
%=========================================================================
function f = logl( p,data )

    f = - sum( loglt( p,data ) );

end
%=========================================================================
%
%   Returns concentrated log-likelihood at each observation
%
%=========================================================================

function f = loglt( p,x )

    y    = x(:,1); 
    ylag = x(:,2);
    st   = x(:,3);
    
    % Transition function
    gam = abs( p(end-1) );
    thr = abs(p(end));
    tmp = (st - thr)/std(st);
    Gt  = 1./( 1+exp( -gam*tmp ) );
 
    % Compute errors
    tmp1 = p(1)+p(2)*ylag; 
    tmp2 = p(3)+p(4)*ylag;
    ut   = y - (tmp1+tmp2.*Gt);
    
    % Concentrate sigma
    s2 = std( ut  );     
   
    % Likelihood function
    f = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*ut.^2/s2;
    
end


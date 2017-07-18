%==========================================================================
%
%   	Program to estimate a bivariate svar model and demonstrate
%   	alternative restrictions
%
%==========================================================================

function svar_bivariate( )
    
    clear all
    clc
    
    % Load data
        % Data from 1970Q1 to 2002Q2    
		%1. Oil price 
		%2. Output EMU
		%3. CPI EMU
		%4. Interest rate EMU
		%5. Output US
		%6. CPI US
        %7. Interest rate US


    load peersman;
    
    % Choose identification method
    % 0 = non-orthogonal shocks
    % 1 = orthogonal shocks with short-run restrictions
    % 2 = orthogonal shocks with long-run restrictions
    % 3 = orthogonal shocks with short- and long-run restrictions
    % 4 = orthogonal shocks with sign restrictions
    itype = 4;      

    % Data transformations
    loil = log( ytdata(:,1) );     %#ok<NASGU,NODEF> % log price of oil    
    lo   = log( ytdata(:,5) );     % log output          
    lp   = log( ytdata(:,6) );     % log price level     
    r    = ytdata(:,7);            %#ok<NASGU> % interest rate       

    % Construct variables for use in VAR: 
    % % growth rates in oil, output and price, and level of interest rate
 
    yact = [ lo lp ];
    yact = trimr(yact,36,0);       % Sample period 1970Q1 to 2002Q2

    y = 400*( trimr(yact,1,0) - trimr(yact,0,1) ) ;

    % Set SVAR parameters
    p = 2;      % Order of VAR lags
    q = 40;     % Order of VMA lags      


    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1) ;
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);

    disp( ' Covariance matrix of VAR residuals ');
    disp( omegav );


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

    % Identification
    if itype == 0       % Non-orthogonal one unit impulse responses 
        s = [ 1 0 ; 0 1 ];  
    elseif itype == 1   % Orthogonal one-standard deviation impulse response based on short-run restrictions
        s = [ 3 0 ; -1 2 ];   
    elseif itype == 2   % Orthogonal one-standard deviation impulse response based on long-run restrictions
        f = [ 3 0  ;
             -1 2 ];
        s = a1*f;          
    elseif itype == 3  % Orthogonal impulse response based on short-run and long-run restrictions
        s11 = 2; 
        s12 = 1; 
        s21 = 0; 
        s22 = (-a1inv(1,1)/a1inv(1,2))*s12; 
        s   = [s11  s12 ; 
               s21  s22 ];   
    elseif itype == 4   % Orthogonal impulse response based on sign short-run restrictions 
        s = [ 3 0;  -1 2];                
        th = pi*0.4;
        qmat = [ cos(th)  -sin(th) ;
                 sin(th)   cos(th) ];
    
        sq = s*qmat';
    
        disp( '------------------------')
        disp( 's' )
        disp( s )
        disp( '------------------------')
        disp( 'qmat' )
        disp( qmat )
        disp( '------------------------')
        disp( 'sq')
        disp( sq )
        s = sq;
    
    end
    
    disp( '------------------------')
    disp( 'a1' )
    disp( a1 )
    disp( '------------------------')
    disp( 'a1inv' )
    disp( a1inv )
    disp( '------------------------')
    disp( 's' )
    disp( s )
    disp( '------------------------')

    
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
    impulse = vec(s');

    for i=2:q+1

       tmp      = ( reshapeg( c(:,i),k,k )*s )';
       impulse  = [ impulse tmp(:) ];

    end 
	impulse = impulse';
 
    % Compute cumulative sum of impulse responses for levels 
    impulse = cumsum( impulse );     

    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    t = 0:1:q;

    %--------------------------------------------------------%
    % Panel (a)
    subplot(2,2,1)
    plot(t,[impulse(1:q+1,1) zeros(q+1,1)],'-k','LineWidth',1);
    if itype == 0
        title( '$v_{1,t}$ shock '); 
    else
        title( '$z_{1,t}$ shock '); 
    end
    ylabel( '$y_{1,t}$' );
    xlabel('Quarter');
    box off
    %axis tight
    
    %--------------------------------------------------------%
    % Panel (b)
    subplot(2,2,2)
    plot(t,[impulse(1:q+1,2) zeros(q+1,1)],'-k','LineWidth',1);
    if itype == 0
        title( '$v_{2,t}$ shock '); 
    else
        title( '$z_{2,t}$ shock '); 
    end
    %ylabel( '$y_{1,t}$' );
    xlabel('Quarter');
    box off
    %axis tight


    %--------------------------------------------------------%
    % Panel (c)
    subplot(2,2,3)
    plot(t,[impulse(1:q+1,3) zeros(q+1,1)],'-k','LineWidth',1);
%     if itype == 0
%         title( '$v_{1,t}$ shock '); 
%     else
%         title( '$z_{1,t}$ shock '); 
%     end
    ylabel( '$y_{2,t}$' );
    xlabel('Quarter');
    box off
    %axis tight


    %--------------------------------------------------------%
    % Panel (d)
    subplot(2,2,4)
    plot(t,[impulse(1:q+1,4) zeros(q+1,1)],'-k','LineWidth',1);
%     if itype == 0
%         title( '$v_{2,t}$ shock '); 
%     else
%         title( '$z_{2,t}$ shock '); 
%     end
%    ylabel( '$y_{2,t}$' );
    xlabel('Quarter');
    box off
    %axis tight



    % Print the tex file to the relevant directory
    %laprint(1,'svar_biv0','options','factory');
    %laprint(1,'svar_biv1','options','factory');
    %laprint(1,'svar_biv2','options','factory');
    %laprint(1,'svar_biv3','options','factory');
    %laprint(1,'svar_biv4','options','factory');

end
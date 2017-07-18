%==========================================================================
%
%   	Program to estimate the Peerman SVAR model
%
%==========================================================================

function svar_peersman( )
    
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

    % Data transformations
    loil = log( ytdata(:,1) );     %#ok<NODEF> % log price of oil    
    lo   = log( ytdata(:,5) );     % log output          
    lp   = log( ytdata(:,6) );     % log price level     
    r    = ytdata(:,7);            % interest rate       

    % Construct variables for use in VAR: 
    % % growth rates in oil, output and price, and level of interest rate
 
    yact = [ loil lo lp r ];
    yact = trimr(yact,36,0);       % Sample period 1970Q1 to 2002Q2

    y = [100*(trimr(yact(:,1:3),1,0) - trimr(yact(:,1:3),0,1))  trimr(yact(:,4),1,0) ];

    % Set SVAR parameters

    p = 3;      % Order of VAR lags
    q = 40;     % Order of VMA lags      


    % Estimate the VAR with p lags, a constant and a time trend     
    ylag   = [ ones(length(y),1) (1:1:length(y))' ];
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);



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
            12.9662120731029020 
            0.0133377541945983 
            0.0917703235318815 
            0.1492472650864925 
            0.3093300004256599 
            -0.1277889658083315 
            -0.1312058376209482 
            0.4307222392355360 
            0.1183146060546325 
            0.6944925202098823
        ];

  options = optimset(   'Display',              'off', ...
                        'MaxIter',              2000,   ...
                        'MaxFunEvals',          4000 ); %,   ...
                        %'TolX',                 1e-8,   ...
                        %'TolFun',               1e-8   ...
                                      
        
            
    [ b,fval ]  = fminunc( @(b) neglogl(b,v,a1inv),bstart,options );     
  
    disp(' Log-likelihood function ');
    disp( fval );

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
    disp(' S matrix at optimal parameters');
    s43 = -(a1inv(2,2)*b(8)+a1inv(2,3)*b(9))/a1inv(2,4);  
    s34 = -(a1inv(2,4)*b(10))/a1inv(2,3);     

    s   =     [   b(1)   0       0       0       ;
        	      b(2)   b(5)    b(8)    0       ;
	              b(3)   b(6)    b(9)    s34     ;
   		          b(4)   b(7)    s43     b(10) ] ;
              
    disp( s )

    impulse = vec(s');

    for i=2:q+1

       tmp      = ( reshapeg( c(:,i),k,k )*s )';
       impulse  = [ impulse tmp(:) ];

    end 
	impulse = impulse';
 
    % Compute cumulative sum of impulse responses for levels ie first three series **/
    impulse = [ cumsum(impulse(:,1:12)) impulse(:,13:16) ];     


    % Construct variance decompositions 
    tmp0   = reshapeg( cumsum( impulse(1:end-1,:).^2 ),q*k,k );
    tmp1   = repmat( sum(tmp0,2),1,k );
    decomp = reshapeg( 100*tmp0 ./ tmp1 , q , k^2 );
 

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
    subplot(4,4,1)
    plot(t,[impulse(1:q+1,1) zeros(q+1,1)],'-k','LineWidth',1);
    title('Oil Price Shock');
    ylabel('Oil price');
    xlabel('Quarter');
    box off
    axis tight
    
    %--------------------------------------------------------%
    % Panel (b)
    subplot(4,4,2)
    plot(t,[impulse(1:q+1,2) zeros(q+1,1)],'-k','LineWidth',1);
    title('Supply Shock');
    %ylabel('Oil price');
    xlabel('Quarter');
    box off
    axis tight


    %--------------------------------------------------------%
    % Panel (c)
    subplot(4,4,3)
    plot(t,[impulse(1:q+1,3) zeros(q+1,1)],'-k','LineWidth',1);
    title('Demand Shock');
    %ylabel('Oil price');
    xlabel('Quarter');
    box off
    axis tight


    %--------------------------------------------------------%
    % Panel (d)
    subplot(4,4,4)
    plot(t,[impulse(1:q+1,4) zeros(q+1,1)],'-k','LineWidth',1);
    title('Money Shock');
    %ylabel('Oil price');
    xlabel('Quarter');
    box off
    axis tight


    %--------------------------------------------------------%
    % Panel (e)
    subplot(4,4,5)
    plot(t,[impulse(1:q+1,5) zeros(q+1,1)],'-k','LineWidth',1);
    ylabel('Output');
    xlabel('Quarter');
    box off
    axis tight


    %--------------------------------------------------------%
    % Panel (f)
    subplot(4,4,6)
    plot(t,[impulse(1:q+1,6) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Output');
    xlabel('Quarter');
    box off
    axis tight

    %--------------------------------------------------------%
    % Panel (g)
    subplot(4,4,7)
    plot(t,[impulse(1:q+1,7) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Output');
    xlabel('Quarter');
     box off
    axis tight

    %--------------------------------------------------------%
    % Panel (h)
    subplot(4,4,8)
    plot(t,[impulse(1:q+1,8) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Output');
    xlabel('Quarter');
    box off
    axis tight

    %--------------------------------------------------------%
    % Panel (i)
    subplot(4,4,9)
    plot(t,[impulse(1:q+1,9) zeros(q+1,1)],'-k','LineWidth',1);
    ylabel('Price');
    xlabel('Quarter');
    box off
    axis tight
    
    %--------------------------------------------------------%
    % Panel (j)
    subplot(4,4,10)
    plot(t,[impulse(1:q+1,10) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Price');
    xlabel('Quarter');
    box off
    axis tight
    
    %--------------------------------------------------------%
    % Panel (k)
    subplot(4,4,11)
    plot(t,[impulse(1:q+1,11) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Price');
    xlabel('Quarter');
    box off
    axis tight
    
    %--------------------------------------------------------%
    % Panel (l)
    subplot(4,4,12)
    plot(t,[impulse(1:q+1,12) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Price');
    xlabel('Quarter');
    box off
    axis tight
    
    
    %--------------------------------------------------------%
    % Panel (m)
    subplot(4,4,13)
    plot(t,[impulse(1:q+1,13) zeros(q+1,1)],'-k','LineWidth',1);
    ylabel('Interest rate');
    xlabel('Quarter');
    box off
    axis tight
     
    %--------------------------------------------------------%
    % Panel (n)
    subplot(4,4,14)
    plot(t,[impulse(1:q+1,14) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Interest rate');
    xlabel('Quarter');
    box off
    axis tight

    %--------------------------------------------------------%
    % Panel (o)
    subplot(4,4,15)
    plot(t,[impulse(1:q+1,15) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Interest rate');
    xlabel('Quarter');
    box off
    axis tight
 
    %--------------------------------------------------------%
    % Panel (p)
    subplot(4,4,16)
    plot(t,[impulse(1:q+1,16) zeros(q+1,1)],'-k','LineWidth',1);
    %ylabel('Interest rate');
    xlabel('Quarter');
    box off
    axis tight

    
    % Print the tex file to the relevant directory
    %laprint(1,'peersman','options','factory');



end
%
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Wrapper function for log-liklihood
%--------------------------------------------------------------------------

function lf = neglogl( b,v,a1 )

    lf = - sum( loglt( b,v,a1 ) );
end

%--------------------------------------------------------------------------
% Log-likelihood function for SVAR
%--------------------------------------------------------------------------

function lf = loglt( b,v,a1inv )

    [ t,n ] = size( v );
    lf      = zeros( t,1 );
    
    % Long-run restrictions
    s43 = -(a1inv(2,2)*b(8)+a1inv(2,3)*b(9))/a1inv(2,4);  
    s34 = -(a1inv(2,4)*b(10))/a1inv(2,3);     

    s   =     [   b(1)   0       0       0       ;
        	      b(2)   b(5)    b(8)    0       ;
	              b(3)   b(6)    b(9)    s34     ;
   		          b(4)   b(7)    s43     b(10) ] ;
			  
	omegav = s*s';

    for i = 1:t;

        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(i,:)*inv(omegav)*v(i,:)';    
    end
    
    
end
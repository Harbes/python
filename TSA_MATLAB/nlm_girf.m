%=========================================================================
%
%     Program to compute generalized impulse response function 
%     (Koop, Pesaran and Potter (1996), Journal of Econometrics)
%
%=========================================================================
function nlm_girf( )

    clear all
    clc
    
    % Initialise the random number generator  
    RandStream.setGlobalStream( RandStream('mt19937ar','seed',12) ); 	

    % Parameters
    delta = -1;             % Size of the shock     
    rho1  = 0.25;           % AR(1) parameter                  
    rho2  = 0.50;           % TAR parameter                    
    rho   = rho1 + rho2;
    t     = 1000;           % Choose sample size       
    n     = 10;             % Forecast horizon of impulses     

    % Simulate the data   
    v = randn(t,1);
    y = zeros(t,1);

    for i = 2:t

        y(i) = rho1*y(i-1) + rho2*y(i-1)*(y(i-1)>=0) + v(i);

    end
    
    % Compute generalized impulse response function 
    % average over horizon within loop for a particular history  
    % and then average over histories
    
    % Total number of draws needed for a given initial condition: 
    % t-1 from history and n+1 from number of impulse horizons       
    tn = (t-1)*(n+1);                                          
    
    impulse = zeros(t-1,n+1);

    % Loop through the data to change the initial condition (history)   
    for i=1:t-1
        
        %  Bootstrap residuals 
        ind    = fix( rand(tn,1)*(t-1) + 1 );
        v_boot = v(ind);              
        v_boot = reshape(v_boot,n+1,t-1);

        ye0 = zeros(n+1,t-1);
        ye1 = zeros(n+1,t-1);
    
        % Loop through horizon of impulse
        for j = 1:t-1                                                                

            % Initial condition based on a boostrap draw              
            ye0(:,j) = model(v_boot(:,j),v_boot(1,j),rho1,rho2);        
            
            % Initial condition based on history (i subscript) plus 1  
            ye1(:,j) = model(v_boot(:,j),v(i)+delta,rho1,rho2);           

        end

        % Average over horizon given an initial condition (history)          
        impulse(i,:) = (mean(ye1,2) - mean(ye0,2))';          
    end
    
    % Average each impulse across histories (ie initial conditions)     
    impulse_girf = mean(impulse)';    
    
    % Linear impulse response    
    impulse_linear = delta*recserar( zeros(11,1) , 1 , rho );   

    format short
	disp('     GIRF     Linear ')
    disp( [impulse_girf impulse_linear] );
    
end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Model used to compute GIRF
%-------------------------------------------------------------------------
function y = model(v,y0,rho1,rho2)

 
    y = y0 + zeros(length(v),1);

    for i=2:length(y)

        % Linear model 
%         rho = rho1+rho2;
%         y(i) = rho*y(i-1) + v(i);           

        % Bilinear model     
        y(i) = rho1*y(i-1) + rho2*y(i-1)*(y(i-1)>=0) + v(i);

    end
end


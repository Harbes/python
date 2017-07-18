%=========================================================================
%
%   Compute GMM estimates of the parameters of a student t distribution
%
%=========================================================================
function gmm_student( )

    clear all;
    clc;
    format short;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) ); 
    
    t = 10;                             
    %y = round(5 + 2*randn(t,1));     
    
    y = load('table.dat','-ascii'); 


    % Using 2 moment conditions 
    % Zero iteration
    theta0 = [mean(y) ; 5];                
    g = numgrad(@gmmcrit2,theta0,y)';           
    h = numhess(@gmmcrit2,theta0,y);

    disp('Zero iteration');
    disp([' Parameter estimate = ', num2str(theta0') ] );
    disp([' Objective function = ', num2str( gmmcrit2(theta0,y) ) ] );
    disp([' Gradient           = ', num2str(g') ] );
    disp(' Hessian            = ');
    disp( h );
    % Newton Raphson update
    theta1 = theta0 - inv(h)*g;     

    
    % First iteration
    g = numgrad(@gmmcrit2,theta1,y)';           
    h = numhess(@gmmcrit2,theta1,y);
 
    disp(' ');
    disp('First iteration');
    disp([' Parameter estimate = ', num2str(theta1') ] );
    disp([' Objective function = ', num2str( gmmcrit2(theta1,y) ) ] );
    disp([' Gradient           = ', num2str(g') ] );
    disp(' Hessian            = ');
    disp( h );
    % Newton Raphson update
    theta2 = theta1 - inv(h)*g;     

    % Second iteration
    g = numgrad(@gmmcrit2,theta2,y)';           
    h = numhess(@gmmcrit2,theta2,y);
 
    disp(' ');
    disp('Second iteration');
    disp([' Parameter estimate = ', num2str(theta2') ] );
    disp([' Objective function = ', num2str( gmmcrit2(theta2,y) ) ] );
    disp([' Gradient           = ', num2str(g') ] );
    disp(' Hessian            = ');
    disp( h );
    
    v = inv(h)/t;
    disp('Covariance matrix   = ');
    disp( v );

    % Iterative solution
    ops = optimset('LargeScale','off','Display','off');
    [thetahat,fc,~,~,~,H]  = fminunc(@(theta) gmmcrit2(theta,y),theta0,ops);
   
    disp(' ');
    disp('Two moment conditions')
    disp([' Parameter estimate = ', num2str(thetahat') ] );
    disp([' Objective function = ', num2str( fc ) ] );
    v = inv(H)/t;
    disp('Covariance matrix   = ');
    disp( v );
 
 
    % Using 3 moment conditions
    [thetahat,fc,~,~,~,H]  = fminunc(@(theta) gmmcrit3(theta,y),theta0,ops);

    disp(' ');
    disp('Three moment conditions')
    disp([' Parameter estimate = ', num2str(thetahat') ] );
    disp([' Objective function = ', num2str( fc ) ] );
    v = inv(H)/t;
    disp('Covariance matrix   = ');
    disp( v ); 

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% GMM objective function - 2 moment conditions
%-------------------------------------------------------------------------
function q = gmmcrit2(theta,y)


        mu = theta(1);
        nu = theta(2);
        m1 = y - mu;
        m2 = (y - mu).^2 - nu/(nu - 2);
        m  = [m1  m2];
        w  = m'*m/length(y);
        q  = 0.5*mean(m)*inv(w)*mean(m)';
end
%-------------------------------------------------------------------------
% GMM objective function - 3 moment conditions
%-------------------------------------------------------------------------
function q = gmmcrit3(theta,y)

        mu = theta(1);
        nu = theta(2);
        m1 = y - mu;
        m2 = (y - mu).^2 - nu/(nu - 2);
        m3 = (y - mu).^4 - 3*nu^2/((nu - 2)*(nu - 4));
        m  = [m1  m2 m3];
        w  = m'*m/length(y);
        q  = 0.5*mean(m)*inv(w)*mean(m)';
end

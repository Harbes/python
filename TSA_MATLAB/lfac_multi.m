%=========================================================================
%
%   Recursions of the multivariate Kalman filter
%
%=========================================================================
function lfac_multi()

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234567) );

    t   = 5;
    % Data from text
    y1 = [1.140, 2.315, -0.054, -1.545, -0.576]; 
    y2 = [3.235, 0.552, -0.689,  1.382,  0.718]; 
    y3 = [1.748, 1.472, -1.413, -0.199,  1.481];
    
    % Data for exercise     
%     y1 = [2.500, 2.017, -0.107, -0.739, -0.992];
%     y2 = [2.000, 1.032, -0.535,  0.061,  0.459]; 
%     y3 = [1.500, 0.047, -0.964,  0.862,  1.910];
    
    y  = [y1'  y2'  y3'];
       
    % Reproduce the simulated GAUSS results 
%    y =   [    2.658    3.490    2.221 ;
%              1.012    0.577   -0.127 ;
%             -1.590   -0.595    2.340 ;
%             -4.941   -4.562   -2.403 ;
%             -4.418   -2.807   -4.605 ];
%         
    % Parameter matrices
    Phi = [ 0.8,  0.0 ;
            0.0,  0.5 ];
    Lam = [ 1.0 ,  0.5 ;
            1.0 ,  0.0 ;
            1.0 , -0.5  ];
    R = [ 0.25  0.00  0.00 ; 
          0.00  0.16  0.00 ; 
          0.00  0.00  0.09];
    Q = eye(2);

    % Kalman filter     
    lnlt(y,Phi,Lam,R,Q);

end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Multivariate Kalman filter
%--------------------------------------------------------------------------
function lnlt(y,Phi,Lam,R,Q)


    % Allocate arrays
    [ t,n ]   = size(y);
    k         = size(Q,1);
    lnl       = zeros(t,1);
    
	% Recursions of the Kalman Filter
    % Initialisation following Harvey and Hamilton  
    st = zeros(k,1);
    pt = reshape(inv(eye(k^2) - kron(Phi,Phi))*Q(:),k,k)';   

    mt = Lam*st;
    vt = Lam*pt*Lam' + R;
    ut = y(1,:)' - mt;

    lnl(1) = - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut'*inv(vt)*ut;
    
    Kal = pt*Lam'*inv(vt);
    s0 = st + Kal*ut;
    p0 = pt - Kal*Lam*pt;

    % Main loop over observations

    for i = 2:t
        
	    % Prediction 
        st = Phi*s0;                     
        pt = Phi*p0*Phi' + Q;	        
              
        % Observation
        mt = Lam*st;
        vt = Lam*pt*Lam' + R;
        ut = y(i,:)' - mt;
       
        % Construct log-likelihood function
        lnl(i) = - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut'*inv(vt)*ut;
                
		% Update    			
        Kal = pt*Lam'*inv(vt);
        s0 = st + Kal*ut;
        p0 = pt - Kal*Lam*pt;
        
        disp('Preditions, Observations and Updates');
        disp('st = ');
        disp(st);
        disp('pt = ');
        disp(pt);
        disp('vt = ');
        disp(vt);
        disp('ut = ');
        disp(ut);
    
        disp('s0 = ');
        disp(s0);
        disp('p0 = ');
        disp(p0);        
    end
end

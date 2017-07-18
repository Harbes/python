%========================================================================
%
%   Estimate the Sims 6-variate SVAR model with long-run restrictions
%
%========================================================================

clear all
clc

% Load data from 1959:1 to 1998:12
load sims_data.mat;      
    
% Define varaibles
r    = ytdata(:,1);        
lex  = log( ytdata(:,2) );
lcp  = log( ytdata(:,3) );
lm   = log( ytdata(:,4) );
lp   = log( ytdata(:,5) );
lo   = log( ytdata(:,6) );
sdum = ytdata(:,7:17);
y    = [r , lex , lcp , lm , lp , lo];   

% arunknown = vgxset('n',6,'nAR',2,'Constant',true);
% [EstSpec,EstStdErrors] = vgxvarx(arunknown,y);
% vgxdisp(EstSpec, EstStdErrors);

% Estimate the VAR with p lags	
p = 14;		
ylag = lagmatrix(y,1:p);
ylag(any(isnan(ylag),2),:) = [];

xmat = [ ones(length(ylag),1)  trimr(sdum,0,p) ylag];
ymat = trimr(y,p,0);

% for i = 1:p  
%     
%     ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];
% 
% end

bar = xmat\ymat;
v   = ymat - xmat*bar;
vc  = v'*v/length(v);

%		Reduced form approach	
u1 = v(:,1);

a2 = v(:,1)\v(:,2);
u2 = v(:,2) - v(:,1)*a2;

a3 = v(:,1:2)\v(:,3);
u3 = v(:,3) - v(:,1:2)*a3;

a4 = v(:,1:3)\v(:,4);
u4 = v(:,4) - v(:,1:3)*a4;

a5 = v(:,1:4)\v(:,5);
u5 = v(:,5) - v(:,1:4)*a5;

a6 = v(:,1:5)\v(:,6);
u6 = v(:,6) - v(:,1:5)*a6;

b0_rf = [   1   ,    0     ,    0    ,    0    ,    0    ,   0  ;
        -a2(1)  ,    1     ,    0    ,    0    ,    0    ,   0  ;
        -a3(1)  , -a3(2)   ,    1    ,    0    ,    0    ,   0  ;
        -a4(1)  , -a4(2)   , -a4(3)  ,    1    ,    0    ,   0  ;
        -a5(1)  , -a5(2)   , -a5(3)  , -a5(4)  ,    0    ,   0  ;
        -a6(1)  , -a6(2)   , -a6(3)  , -a6(4)  , -a6(5)  ,   1  ] ;
    
u     = [u1 , u2 , u3 , u4 , u5 , u6];
d_rf  = u'*u/size(u,1);

% Choleski decomposition approach	
s      = chol(vc)';
b0inv  = bsxfun(@rdivide,s,diag(s)');
b0_cd = inv(b0inv);
d_cd = diagrv(eye(6),diag(s).^2); 

disp('B0: Reduced form')
disp( b0_rf )

disp('B0: Choleski decomposition')
disp( b0_cd )




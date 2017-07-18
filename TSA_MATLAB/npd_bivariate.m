% ========================================================================
%
%    Bivariate kernel  
%
% ========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

% Generate data
n  = 2;                      % Dimension                          
t  = 20000;                   % Sample size
yt = normrnd( 0,1,[t n]);


% Set up a grid between -4 and 4
ymesh = (-4:0.2:4)';

% Need to work out all possible points at which to compute bivariate
% density on this mesh and express results as vectors
[y1m,y2m] = meshgrid(ymesh);

y = [reshape( y1m,[],1 ) reshape( y2m,[],1 )];

% Estimate density using product kernel
fac = 1/(4.0+n);
h   = std( yt )./(t.^fac);
ph  = prod( h );

ker  = zeros(n,1);
pdf  = zeros(length(y),1);
pker = zeros(t,1);
for j = 1:length(y)
    
    for i = 1:t;
        
        for p = 1:n
            ker(p) = normpdf( (y(j,p) - yt(i,p))/h(p) );
        end 
        pker(i) = prod( ker );
     end
     pdf(j) = mean( pker )/ph;
end

%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

mesh(y1m,y2m,reshape( pdf,length(y1m),length(y2m) ),'EdgeColor','black')

xlabel('$y_{1,t}$');
ylabel('$y_{2,t}$');
zlabel('$f(y_{1,t},y_{2,t})$');
set(gca,'ztick',[])
axis tight
grid 'off'
box 'off'


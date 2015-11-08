% Normal Inverse Gamma test

% generate observations
N=500

% The coefficients are [2, 1], with the first one denoting the y intersect
beta=[5; -1];
printf("Create observations on line with coefficients: y=%f+%f*x\n", beta(1), beta(2));

X0=[ones(1,N); sort(unifrnd(0,10,size(beta,1)-1,N))]';
y0=(X0*beta);

%beta=[1; 2];
%X1=[ones(1,N); sort(unifrnd(0,10,size(beta,1)-1,N))]';
%y1=(X1*beta);

%y=[y0; y1];
%X=[X0; X1];
y=y0;
X=X0;

sigma=0.1;

y += normrnd(0, sigma, length(y), 1);

plot(X(:,2), y,'.')

% prior
S0.mu=[1; 1];
S0.Lambda=[0.5 0; 0.5 0];
S0.a=1;
S0.b=2;

z = [X y]';

Sn = nigupdate(z, S0);

% Inferred line
printf("The inferred line has coefficients: y=%f+%f*x\n", Sn.mu(1), Sn.mu(2));


data_dir='../data';

fname='twolines';
%fname='oneline';

output_file=[data_dir '/' fname '.pnts.data.' date];
% generate observations
N=100

r0=-10;
r1=10;

% The coefficients are [2, 1], with the first one denoting the y intersect
beta=[5; -1];
printf("Create observations on line with coefficients: y=%f+%f*x\n", beta(1), beta(2));

X0=[ones(1,N); sort(unifrnd(r0,r1,size(beta,1)-1,N))]';
y0=(X0*beta);

beta=[1; 2];
printf("Create observations on line with coefficients: y=%f+%f*x\n", beta(1), beta(2));
X1=[ones(1,N); sort(unifrnd(r0,r1,size(beta,1)-1,N))]';
y1=(X1*beta);

y=[y0; y1];
X=[X0; X1];

% add class labels
C=[ones(length(y0),1); ones(length(y1),1)*2];

sigma=0.1;

y += normrnd(0, sigma, length(y), 1);

dlmwrite(output_file, [X(:,2), y, C], 'delimiter', '\t', 'precision', 10);

%plot(X(:,2), y,'.')


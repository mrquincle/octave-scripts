% generate samples on a line and get endpoints

lb=-4;
hb=5;
sample_cnt=100;
sample_sub_cnt = [1 3 10 sample_cnt];
sigma=0.5;
avg=4.5;

lim = 10;
axis_limits=[-lim, +lim, -lim, +lim];

X=unifrnd(lb, hb, 1, sample_cnt);

Y=normrnd(avg, sigma, 1, sample_cnt);

P=[X Y];


alpha=10;

endpoint_cnt=10;

fh=figure(1);

% priors, these are now sampled from a Normal distribution
% if they are small, they are noninformative, if they are large they will
% dominate the data (the data will never be able to get the segment smaller)
pA = -normrnd(0,5);
pB = normrnd(0,5);

plot(X, Y, '.');
axis(axis_limits);

% we can use a flat prior
% f(u|y_s) \sim exp^(-1/2sigma^2/n)*(u-y_s)^2

% or we can use a normal prior density
% f(u|y_s) \sim exp^(-1/2sigma^2)*(u-y_s)^2 * exp^(-1/2s^2)*(u-m)^2
% this is for single observation
% we can update m and s instead of calculating it all the time
% posterior precision: 1/s'^2 = (sigma^2 + n*s^2) / (sigma^2 * s^2)
% posterior mean: m' = m/s^2 (n/sigma^2) + ....

N=sample_cnt;

% let us first assume variance known (sigma is std.dev, so sigma^2)

prior_variance = 0.4;
prior_precision = 1/(prior_variance^2);

prior_mean = 1
%empirical_mean = mean(abs(Y-mean(Y)));
empirical_mean=mean(Y);
empirical_variance = mean((Y-mean(Y)).^2);
empirical_precision = 1/empirical_variance;

% https://noppa.aalto.fi/noppa/kurssi/s-114.1310/luennot/luento_12.pdf
% http://jackman.stanford.edu/classes/BASS/ch2.pdf
% http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf

% use normal conjugacy to find out that N observations with variance sigma^2 and mean x'
% is same as single observation x' with variance sigma^2/N

known_sigma = sigma;
known_precision = 1/(sigma^2);
known_mean = avg;

post_precision = (prior_precision + N*known_precision);

post_mean = (prior_mean * prior_precision + empirical_mean * N * known_precision) / post_precision

post_variance = 1 / post_precision;

% note, we actually note the precision, so we should throw away this posterior precision
% which is incorrect

post_sigma = sqrt(post_variance);

% actually, we would like to use a prior from
% http://jakevdp.github.io/blog/2014/06/14/frequentism-and-bayesianism-4-bayesian-in-python/
% p(a,beta,sigma) \sim 1/sigma (1 + beta^2)^(-3/2)

% exponent exp( -100 ) with large negative values are basically 0. Considering this is a sum of squares and we only
% divide by something like 1/sigma^2, this is certainly not gonna work.
function result = likelihood(X, Y, alpha, beta, sigma)
	Y_model = alpha + beta * X;
	N = size(X, 2);
	result = 1/((2*pi*sigma^2)^(N/2)) * exp ( -1/(2*sigma^2) * sum(Y-Y_model).^2 );
end

function result = prior(alpha, beta, sigma)
	result = 1/sigma * (1 + beta^2)^(-3/2);
end

function result = posterior(X, Y, alpha, beta, sigma)
	result = likelihood(X, Y, alpha, beta, sigma) * prior(alpha, beta, sigma);
end

function result = loglikelihood(X, Y, alpha, beta, sigma)
	Y_model = alpha + beta * X;
	N = size(X, 1);
	% TODO: Check this, sum of log seems redundant
	result = -1/2 * sum(log(2*pi*sigma^2) + ((Y-Y_model).^2)/(sigma^2) );
end

function result = logprior(alpha, beta, sigma)
	result = -log(sigma) + log(1 + beta^2) * (-3/2);
end

function result = logposterior(X, Y, alpha, beta, sigma)
	result = loglikelihood(X, Y, alpha, beta, sigma) + logprior(alpha, beta, sigma);
end

alpha = 4.5;
beta = 0;
sigma = 0.5;
logpost = logposterior(X, Y, alpha, beta, sigma);
logpost
exp(logpost)

function result = test(data,x,y)
	result = sum(log(normpdf(data,x,y)));
end

data=randn(100,1)*2+3;
logmodelprior=@(m)0; %use a flat prior.
loglike=@(m)test(data,m(1),m(2));
minit=[0 1];
m=mcmc(minit,loglike,logmodelprior,[.2 .5],10000);
m(1:100,:)=[]; %crop drift
plotmatrix(m);

return;

%pA = -20;
%pB = 20;

for p=1:4
	subplot(2,2,p);

	ssc=sample_sub_cnt(p);

	S=P(1:ssc);
	Z=zeros(1:ssc, 1);
	[A B] = extreme_pnts(S);

	%S0 = rnd_pareto(B, alpha, ssc, 1, endpoint_cnt);
	A = min(A, pA);
	B = max(B, pB);

	[S0 S1] = rnd_pair_two_sided_pareto(A, B, alpha, ssc, 1, endpoint_cnt);

	plot(S, Z, 'b*');
	hold on;
	plot(S0, 0, 'r.');
	plot(S1, 0, 'g.');
	axis(axis_limits);
	hold off;
end

% but what we actually want to do is to define a prior instead of just
% sampling from a Pareto distribution

% so, let's assume we have a Gaussian for S0
% how do we update S0?

% okay, above, solves that!! :-), yeah!!


% http://math.arizona.edu/~jwatkins/f-transform.pdf
% http://research.microsoft.com/en-us/um/people/minka/papers/minka-uniform.pdf

% now, we need to incorporate the spread of points
% suppose the points are spread across a horizontal line (but not necessary 0)
% we would have a prior across the y-axis
% and we would update it using the samples
% this is independent from above, so, pretty sure, this works

% however, now if it is a line at an angle
% we can rotate all points with the angle, such that it is a horizontal line
% we identify a prior for the rotation
% however, now we have an equation for the line
% does that not conflict with the pareto prior? no... we can do the Pareto after
% mapping all points on the 1D line






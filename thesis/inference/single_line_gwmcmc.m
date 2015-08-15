% inference for line

addpath("..")

do_plot=true;

% Uniformly distribute points X on line segment between [lb hb], in a total number of N. Then translate and rotate
% these points according to Y=alpha+beta*X. Last, but not least add noise with variance equal to sigma.
lb=-6;
hb=6;
N=50;
alpha=2;
beta=-4;
sigma=0.1;

printf("Create line y=%f+%f*x with noise %f from x=%f to x=%f\n", alpha, beta, sigma, lb, hb);

X=unifrnd(lb, hb, 1, N);
Y=alpha+beta*X;
Y=Y+normrnd(0, sigma, 1, N);

P=[X; Y];

% If set, plot the input distribution of points
if (do_plot)
	figure(1);
	subplot(2,1,1);
	plot(X, Y, '.');
end

% model-dimension is 5, we want to get alpha, beta, and sigma (if we do not fix sigma, sampling goes haywire)
% plus we have the two endpoints that limit the line to a line segment
dim=5;

% We do not bother about infering the hyperparameters yet. Of course, they have to get a (vague) prior as well.
xi=0.001;
nu=-0.001;
K=3;

% order probability functions in successive order of complexity
logPfuns={
	@(m)logprior_line(m(1), m(2), m(3));
	@(m)logprior_segment(m(4), m(5), xi, nu, K);
	@(m)loglikelihood_line(P, m(1), m(2), m(3));
	@(m)loglikelihood_segment(P, m(4), m(5), m(2));
};

%make a set of starting points for the entire ensemble of walkers
walkers=dim*4;
minit=randn(dim,walkers);

% The MCMC hammer is not allowed to start in a probability region with probability 0!

% Hence, we use several pieces of information to initialize our samplers!!
% 1.) We first make all variables larger than 0. This is for example necessary for sigma
% 2.) We set the positive endpoint far from where we expect to be the positive endpoint of the line
% 3.) We set the negative endpoint far from where we expect to be the negative endpoint of the line
minit=abs(minit);
minit(4,:) = minit(4,:) * 10 + 20;
minit(5,:) = minit(5,:) * -10 - 20;

%Apply the MCMC hammer
[models,logP]=gwmcmc(minit,logPfuns,100000);
models(:,:,1:floor(end/5))=[]; %remove 20% as burn-in
models=models(:,:)'; %reshape matrix to collapse the ensemble member dimension

function result = logTotal(a, b, s, x0, x1, P, xi, nu, K)
	result = logprior_line(a, b, s) + logprior_segment(x0, x1, xi, nu, K) + loglikelihood_line(P, a, b, s) + loglikelihood_segment(P, x0, x1, b);
end

if (do_plot)
	subplot(2,1,2);

	%plot3(models(:,1),models(:,2),models(:,3), '.')
	%plot(models(:,1),models(:,2), 'r.')

	N=size(models,1)
	skip=0;
	done=0;
	% show a few of the MCMC steps
	for m = 1:30:N
		x0 = models(m,4);
		x1 = models(m,5);
		a = models(m,1);
		b = models(m,2);
		s = models(m,3);

		if (logTotal(a,b,s,x0,x1,P,xi,nu,K) < -100)
			skip++;
		else
			plot_line_segment(x0, x1, a, b);
			done++;
		end
		if (done > 20)
			break;
		end
	end
	printf("Skipped: %i\n", skip);
	plot(X, Y, 'b.');
end

Q=prctile(models,[0 5 50 95])

for i = 1:4
	R = logTotal(Q(i,1), Q(i,2), Q(i,3), Q(i,4), Q(i,5), P, xi, nu, K)
end



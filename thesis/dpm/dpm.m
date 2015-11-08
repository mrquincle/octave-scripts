% Dirichlet Process Mixture

% The code originates from Francois Caron from the University of Oxford.
% It has been adjusted for different applications in computer vision.
% Some of the original code might still look awkward, because of this history.
%
% There are now multiple models implemented with different priors/likelihoods.
%  - the Normal inverse Wishart prior for a Normal likelihood function with
%    unknown mean and unknown covariance matrix
%  - the Normal inverse Gamma (NIG) prior for a Normal likelihood function with
%    unknown mean and variance
%
% The choice which prior is used is through the variable type_prior, which
% currently has the options 'NIG' and 'NIW'.
%
% The prior is defined in hyperG0. The NIG prior is used as a prior for the
% coefficients for a line.
%   hyperG0.prior = 'NIG';
%   hyperG0.mu = [2;0];
%   hyperG0.a = 10;
%   hyperG0.b = 0.1;
%   hyperG0.Lambda = [ 1 0.5; 0.1 1];
% Likewise it is possible to define a hyperG0.prior = 'NIW' prior.
%
% There are different MCMC algorithms implemented. These can be tested through
% setting the variable type_algo. The options tested for all cases are
% 'BMQ' and 'CRP'. It seems 'collapsedCRP' is implemented incorrectly.
%
% To test the classification results, RandIndex and similar metrics are used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

addpath('../performance')
addpath('../inference')
addpath('../inference/prior/nig')
addpath('../inference/prior/niw')
addpath('../inference/likelihood/normal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the prior to use
prior = {'NIQ'; 'NIG'};
type_prior=prior{2}
printf('Use prior ''%s''\n', type_prior);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the dataset

outfname='twolines';
data_dir='../data/lines';
%data_dir='/home/anne/myworkspace/aim/mrquincle/CornerDetectorModule/build_standard';
data_glob=[data_dir '/*.data.txt'];

output_dir='../output';
mkdir(output_dir);

timestamp=strtrim(ctime(time()));
usec=num2str(gmtime(now()).usec)
timestamp( timestamp == " " ) = "_";

output_inference_file=[output_dir '/' outfname '.pnts.' timestamp '.' usec '.data.txt'];

fileList = glob(data_glob);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some other configuration options

% Number of iterations
niter = 100;

% Type of plots (plot at every Gibbs step, or only after each point is updated)
doPlot = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet Process Mixture parameters

% The most important parameter in a DPM is the scale or concentration parameter
% called alpha. With a small alpha we will have only a few clusters. With alpha
% very large we will have many clusters.
alpha = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The inference method
algorithm = {'BMQ'; 'CRP'; 'collapsedCRP'; 'slicesampler'};
type_algo = algorithm{2};
printf('Use algorithm ''%s''\n', type_algo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Persistent, so we can perform post-analysis on our inference
save(output_inference_file, '-append', 'output_inference_file', 'niter',
	'doPlot', 'type_prior', 'alpha', 'type_algo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We perform inference over a list of datasets

for f = 1:length(fileList)

	% Clear all variables when we get a new dataset, except for the ones we
	% define here: the output file name, prior type, number of iterations,
	% plotting variable, alpha, and the algorithm type
	clear -x fileList output_inference_file type_prior f niter doPlot alpha type_algo

	% Load the data from the dataset
	data_file=fileList{f,1}
	pnts=load(data_file);

	% Assume last item is label, so remove for training
	P=pnts(:,1:end-1)';

	% Prepend the X-dimension with a constant 1. This is unique for a
	% transformation to [y-intercept slope] coordinates.
	X0=ones(1,size(P,2));
	P=[X0; P];

	% Plot the raw dataset, the input for the algorithm
	if (doPlot != 0)
		y=P(end,:);
		X=P(2:end-1,:);

		figure('name', 'simulated data')
		plot(X, y, '.')
		xlabel('X')
		ylabel('Y')
		xlim([-1 1]*20);
		ylim([-1 1]*20);
	end

	% Set the parameters of the base distribution G0
	switch(type_prior)
	case 'NIG'
		hyperG0.prior = 'NIG';
		% this mean is used as prior for the coefficients for a line, 0, 0 means a horizontal line through the origin
		% [0,1] means a horizontal line through 1, [1,0] means a diagonal line through the origin (y=x)
		hyperG0.mu = [2;0];
		hyperG0.a = 10; % > 0
		hyperG0.b = 0.1;
		hyperG0.Lambda = [ 1 0.5; 0.1 1];
	case 'NIW'
		hyperG0.prior = 'NIW';
		% using column-vectors
		hyperG0.mu = [0;0];
		hyperG0.kappa = 1;
		hyperG0.nu = 4;
		hyperG0.lambda = eye(2);
	otherwise
		error('Unknown prior ''%s''', type_prior);
	end

	%%
	% # Look at the Matlab function gibbsDPM_algo1.
	% # Run the Gibbs sampler for 20 iterations, with graphical output.
	[c_st, c_est, similarity] = gibbsDPM(P, hyperG0, alpha, niter, type_algo, doPlot);

	% Plot posterior similarity matrix (data ordered with respect to point estimate)
	if (doPlot)
		[~, ind] = sort(c_est);
		figure
		imagesc(double(similarity(ind, ind)))

		% Plot point estimate of the partition
		cmap = colormap;
		figure
		for i=1:size(P, 2)
			plot(P(2,i), P(3,i), 'o', 'color', cmap(mod(10*c_est(i),63)+1,:), 'linewidth', 3);
			hold on
		end
		hold off
		drawnow()
		title('Bayesian point estimate')
	end

	% we should somehow be able to compare the assignments with the ground truth
	R=[c_est pnts(:,end)];
	[AR,RI,MI,HI]=RandIndex(R(:,1), R(:,2));

	% for most likely cluster assignment, subsequently establish lines
%	for i=1:c_est(ind(end))
%		Pl=P(:,c_est==i);
%		Pu=update_SS(Pl,hyperG0);
%		number=size(Pl,2)
%		mu=Pu.mu'
%	end

	% find MAP
	for i=1:size(c_st, 2)
		logjoint(i) = logjoint_dpm(P, c_st(:,i), alpha, hyperG0);
	end
	[~, indmap] = max(logjoint);
	c_map = c_st(:, indmap);
	if (doPlot)
		figure;
		plot(logjoint)
		xlabel('MCMC iteration')
		xlabel('log joint distribution')

		% Plot MAP estimate of the partition
		cmap = colormap;
		figure
		for i=1:size(P, 2)
			plot(P(2,i), P(3,i), 'o', 'color', cmap(mod(10*c_map(i),63)+1,:), 'linewidth', 3);
			hold on
		end
		title('MAP estimate')
	end

	% we should somehow be able to compare the assignments with the ground truth
	RT=[c_map pnts(:,end)];
	[AR,RI,MI,HI]=RandIndex(RT(:,1), RT(:,2))

	clen=size(unique(c_map),1);
	chist=zeros(2, clen);
	[chist(1,:) chist(2,:)] = hist(c_map, unique(c_map));

	mus=zeros(size(hyperG0.mu, 1), clen);
	for i=1:clen
		% subset of points belonging to non-empty cluster chist[i]
		Pl=P(:,c_map==chist(2,i));
		Pu=update_SS(Pl,hyperG0);
		mus(:,i)=Pu.mu';
	end
	save(output_inference_file, '-append', 'data_file', 'chist', 'AR', 'RI', 'MI', 'HI', 'mus', 'RT');
end

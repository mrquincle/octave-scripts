%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to generate N lines with a specified number of outliers and Gaussian noise
%
% Author: Anne van Rossum
% Date: Jul, 2014
% Copyrights: Almende B.V., Rotterdam, The Netherlands
% License: LGPLv3
%
% To create these scripts, the following teaching material from Caron has been used:
% http://www.stats.ox.ac.uk/~caron/code/abs2014/
% http://www.stats.ox.ac.uk/~caron/code/abs2014/html/BNP_clustering.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load configuration from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

do_plot=false;
display_avg_cov=false;

fname='dpmline4';
config_dir='config';

data_dir='data/lines';
% Create data directory
mkdir(data_dir);

input_file=[config_dir '/' fname '_config.m'];
run(input_file);

% Number of files to generate
nfiles=100;
for file=1:nfiles

timestamp=strtrim(ctime(time()));
usec=num2str(gmtime(now()).usec)
timestamp( timestamp == " " ) = "_";

output_data_file=[data_dir '/' fname '.pnts.' timestamp '.' usec '.data.txt']
output_config_file=[data_dir '/' fname '.pnts.' timestamp '.' usec '.config.txt']

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default configuration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure Matlab-compatibility by casting broadcast warnings into errors
warning ('error', 'Octave:broadcast');

% include path for the Chinese Restaurant Process algorithm and the Normal-Inverse Wishart Random generator
addpath('/home/anne/myworkspace/octave/nonparam_workshop')
addpath('/home/anne/myworkspace/octave/nonparam_workshop/private')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We run a Chinese Restaurant Process (CRP) for n items resulting a Dirichlet Process (DP) distribution. This is a
% function that calculates the probabilities for all items in a batch. So, we run it only once. The resulting partition
% is fixed.
[partitions clustersize clusters partition_bins] = crprnd(alpha,n);

% This means there is a fixed number of clusters rolling out of the CRP process above.
S=clusters;
fprintf('Number of clusters: %i\n', clusters);

% Initialize the parameters to the right dimensions
%  mu: mean
%  covar: covariance
%  lp: line endpoint
%  lp_covar: not used
%  P: points on the calculated line segment
%  Y: assignments
mu=zeros(dim,S);
covar=zeros(dim,dim,S);
sli=zeros(dim,S);
lp=zeros(dim,S);
lp_covar=zeros(dim,dim,S);
P=zeros(n,dim);
Y=zeros(n,dim);

% We sample for each component from two different latent spaces. One corresponding to the center of a line segment, the
% other corresponding to its angle and length.

% For each cluster we get a new 'mean' vector and 'covariance' matrix for the center of the line segment.
% For each cluster we also get a new 'mean' vector and 'covariance' matrix for each (angle, segment size)
for i=1:S
	% We get a mean and covariance matrix at once through the normal-inverse Wishart random generator
	[mu(:,i) covar(:,:,i)]=normalinvwishrnd(hyper);

	% We get a mean and covariance matrix for the endpoint, but will throw away the covariance matrix later
	% This needs attention! The proper prior for an angle is different!
	[lp(:,i) lp_covar(:,:,i)]=normalinvwishrnd(hyper);

	% Same angle
	if (test_same_angle)
		lp(:,i) = [same_angle ; lp(2,i)];
	end

	% Same length
	if (test_same_length)
		lp(:,i) = [lp(1,i); same_length];
	end
end

% For each pair of (angle,distance) calculate end point
lp(:,:)=[ (lp(2,:).*cos(lp(1,:))); (lp(2,:).*sin(lp(1,:)))];

% Transform to slope-intercept
slope= lp(2,:) ./ lp(1,:);
sli(:,:)=[slope; mu(2,:) - mu(1,:).*slope];

% Display some stuff for debugging purposes, only use if you have a few clusters
if (display_avg_cov)
	fprintf('The averages are:\n')
	mu
	fprintf('And the covariance matrices:\n')
	covar
end

% Here we generate points on the line segment according to a certain distribution. Moreover, we use the covariance
% from the inverse Wishart to account for the spread around that point picked on the line segment. Note that if
% picking a point from a line segment is from a normal distribution, this is not done properly yet. For starters
% points can be beyond the end of the line segment with probability > 0.

% Loop over all data items
for i=1:n
	% Get the cluster id (which is the partition index)
	c=partitions(i);

	% Line distribution
	switch (noise_distribution)
	case 'normal'
		% points on line segment are generated normally
		sn=randn(1);
	case 'uniform'
		% points on line segment are generated uniformly
		sn=2*rand(1)-1;
	end

	% We use a multiplicative structure to pick a point on the line segment
	% The end of the line (lp) compared to the mean is taken negative as well as positive to form a line segment
	P(i,:)=mu(:,c) + lp(:,c)*sn + chol(covar(:,:,c))' * randn(dim,1);
	Y(i,:)=c;
end

P=P';
Y=Y';

if (do_plot)
	switch (dim)
	case 2
		% In the 2D case we plot a bit more, e.g. the 'true' lines that generate the points as well as the mean values.
		plot(P(1,:),P(2,:),'.')
		hold on
		% Print mean values as a red circle
		plot(mu(1,:),mu(2,:),'or');
		% The factor 'f' indicates where exactly the endpoint ends (think of the normal distribution)
		% The factor 'e' indicates where a triangle symbol will be placed.
		% They have nothing to do with the actual data and are just visually helping you to understand
		switch (noise_distribution)
		case 'normal'
			f=1.3;
			g=1;
		case 'uniform'
			f=1; g=1;
		end
		plot(mu(1,:)-g*lp(1,:),mu(2,:)-g*lp(2,:),'^r');
		plot(mu(1,:)+g*lp(1,:),mu(2,:)+g*lp(2,:),'vr');
		for i=1:S
			plot([mu(1,i)-f*lp(1,i); mu(1,i)+f*lp(1,i)],[mu(2,i)-f*lp(2,i), mu(2,i)+f*lp(2,i) ],'-');
		end
		for i=1:S
			scale=20;
			plot([-1,1]*scale,[sli(2,i)+sli(1,i)*-scale, sli(2,i) + sli(1,i)*scale ],'g-');
		end
		xlabel('X');
		ylabel('Y');
		xlim([-1 1]*20);
		ylim([-1 1]*20);
	case 3
		plot3(P(1,:),P(2,:),P(3,:),'.')
		hold on
		% Print mean values as a red circle
		plot3(mu(1,:),mu(2,:),mu(3,:),'or');
	end
end

switch(dim)
case 2
	% prepend file with meta information
	save(output_config_file, 'alpha', 'hyper', 'mu', 'lp', 'sli', 'clusters')
	dlmwrite(output_data_file, [P(1,:)' P(2,:)' Y(1,:)'], 'delimiter', '\t', 'precision', 10);
case 3
	dlmwrite(output_data_file, [P(1,:)' P(2,:)' P(3,:)', Y(1,:)], 'delimiter', '\t', 'precision', 10);
end

% end for loop
end

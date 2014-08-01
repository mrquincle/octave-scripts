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

fname="dpmline2";
config_dir="config";

data_dir="data";
% Create data directory
mkdir(data_dir);

input_file=[config_dir "/" fname ".config"];
output_file=[data_dir "/" fname ".pnts.data"];

source(input_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default configuration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure Matlab-compatibility by casting broadcast warnings into errors
warning ("error", "Octave:broadcast");

% include path for the Chinese Restaurant Process algorithm and the Normal-Inverse Wishart Random generator
addpath("/home/anne/myworkspace/octave/nonparam_workshop")
addpath("/home/anne/myworkspace/octave/nonparam_workshop/private")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We run a Chinese Restaurant Process (CRP) for n items resulting a Dirichlet Process (DP) distribution. This is a 
% function that calculates the probabilities for all items in a batch. So, we run it only once. The resulting partition 
% is fixed.
[partitions clustersize clusters partition_bins] = crprnd(alpha,n);

% This means there is a fixed number of clusters rolling out of the CRP process above.
S=clusters;
printf("Number of clusters: %i\n", clusters);

% Initialize the parameters to the right dimensions
%  mu: mean
%  covar: covariance
%  lp: line endpoint
%  lp_covar: not used
%  P: points on the calculated line segment
mu=zeros(dim,S);
covar=zeros(dim,dim,S);
lp=zeros(dim,S);
lp_covar=zeros(dim,dim,S);
P=zeros(n,dim);

% For each cluster we get a new "mean" vector and "covariance" matrix for the center of the line segment.
% For each cluster we also get a new "mean" vector and "covariance" matrix for one endpoint of the same line segment.
for i=1:S
	% We get a mean and covariance matrix at once through the normal-inverse Wishart random generator
	[mu(:,i) covar(:,:,i)]=normalinvwishrnd(hyper);
	% We get a mean and covariance matrix for the endpoint, but will throw away the covariance matrix later
	[lp(:,i) lp_covar(:,:,i)]=normalinvwishrnd(hyper);

	% Same angle
	%lp(:,i) = [pi ; lp(2,i)];
	
	% Same length
	lp(:,i) = [lp(1,i); 4];
end

% For each pair of (angle,distance) calculate end point
lp(:,:)=[ (lp(2,:).*cos(lp(1,:))); (lp(2,:).*sin(lp(1,:)))];


% Remark: the way the line segments are build up in this way does not lend itself to a prior that for example puts 
% emphasis on vertical or horizontal lines above lines in other directions. For that to exist, it is important to have
% the angle explicitly incorporated in the generative process.

% The use of a Gaussian Process to generate a function is actually to get an independent Gaussian into the game.
% A linear function with points distributed in a Gaussian fashion over it. Is just not the same as generating everything 
% from a 4D Gaussian with dimensions (center (2D), length (1D), angle (1D)). The points on the same line correspond to
% exactly the same center, angle, length parameters and the noise is only in the observation space, not in the 
% parameter space. In the parameter space is it is uncertainty, not noise. The uncertainty is represented by a uniform
% distribution of center, angle, and length (limited by the size of the image). The noise is represented by a normal
% distribution.

% Plot some stuff for debugging purposes, only use if you have a few clusters
if (plot_avg_cov)
	printf("The averages are:\n")
	mu 
	printf("And the covariance matrices:\n")
	covar
end

% Loop over all data items 
for i=1:n
	% Get the cluster id (which is the partition index)
	c=partitions(i);

	% The end of the line (lp) compared to the mean is taken negative as well as positive to form a line segment

	switch (noise_distribution) 
	case "normal"
		% points on line segment are generated normally
		P(i,:)=mu(:,c) + lp(:,c)*randn(1) + chol(covar(:,:,c))' * randn(dim,1);
	case "uniform"
		% points on line segment are generated uniformly
		P(i,:)=mu(:,c) + lp(:,c)*(2*rand(1)-1) + chol(covar(:,:,c))' * randn(dim,1);
	endswitch
end

P=P';

switch (dim)
case 2
	% In the 2D case we plot a bit more, e.g. the "true" lines that generate the points as well as the mean values.
	plot(P(1,:),P(2,:),'.')
	hold on
	% Print mean values as a red circle
	plot(mu(1,:),mu(2,:),'or');
	% The factor "f" indicates where exactly the endpoint ends (think of the normal distribution)
	% The factor "e" indicates where a triangle symbol will be placed. 
	switch (noise_distribution)
	case "normal"
		f=1.3;
		g=1;
	case "uniform"
		f=1; g=1;
	endswitch
	plot(mu(1,:)-g*lp(1,:),mu(2,:)-g*lp(2,:),'^r');
	plot(mu(1,:)+g*lp(1,:),mu(2,:)+g*lp(2,:),'vr');
	for i=1:S
		plot([mu(1,i)-f*lp(1,i); mu(1,i)+f*lp(1,i)],[mu(2,i)-f*lp(2,i), mu(2,i)+f*lp(2,i) ],'-');
	end
	dlmwrite(output_file, [P(1,:)' P(2,:)'], '\t', "precision", 10);
case 3
	plot3(P(1,:),P(2,:),P(3,:),'.')
	hold on
	% Print mean values as a red circle
	plot3(mu(1,:),mu(2,:),mu(3,:),'or');
	dlmwrite(output_file, [P(1,:)' P(2,:)' P(3,:)'], '\t', "precision", 10);
endswitch



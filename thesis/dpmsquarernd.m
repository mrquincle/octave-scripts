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

fname="dpmsquare1";
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
printf("Number of clusters: %i\n", S);

% Initialize the parameters to the right dimensions
%  mu: mean
%  covar: covariance
%  lp: line endpoint
%  lp_covar: not used
%  P: points on the calculated line segment
mu=zeros(dim,S);
covar=zeros(dim,dim,S);

ls=zeros(hyper1.dim,S);

lp1=zeros(dim,S);
lp1_covar=zeros(dim,dim,S);
lp2=zeros(dim,S);
lp2_covar=zeros(dim,dim,S);
P=zeros(n,dim);

% We add an additional CRP to get the same angle for line segments with a probability larger than 0, indicating 
% correlations coming from higher-level objects, such as squares and alike.
[s_p s_cs s_c s_pb] = crprnd(alpha1,S);
SC=s_c;
printf("Number of super clusters: %i\n", SC);

% Normally distributed
s_mu=zeros(hyper1.dim,SC);
s_covar=zeros(hyper1.dim,hyper1.dim,SC);

for i=1:SC
	[s_mu(:,i) s_covar(:,:,i)]=normalinvwishrnd(hyper1);
end

% We sample for each component from two different latent spaces. One corresponding to the center of a line segment, the
% other corresponding to its angle and length.

% For each cluster we get a new "mean" vector and "covariance" matrix for the center of the line segment.
% For each cluster we also get a new "mean" vector and "covariance" matrix for each (angle, segment size)
for i=1:S
	% We get a mean and covariance matrix at once through the normal-inverse Wishart random generator
	[mu(:,i) covar(:,:,i)]=normalinvwishrnd(hyper);

	% We get a mean and covariance matrix for the endpoint, but will throw away the covariance matrix later
	% This needs attention! The proper prior for an angle is different!
	%[lp(:,i) lp_covar(:,:,i)]=normalinvwishrnd(hyper);

	% Get the super-cluster
	sc = s_p(i);
	% Get the angle and two lengths for 2D, get three angles and two lengths for 3D
	ls(:,i)=s_mu(:,sc);
end

% For each pair of (angle,distance) calculate end point
switch(dim) 
case 2
	% We assume order (angle, length1, length2), we make a "square" from it, actually a parallelogram by picking two
	% other points. Then all points will be between the middle \mu and points \mu+-lp1 and \mu+-lp2.
	lp1(:,:)=[ 	(ls(2,:).*cos(ls(1,:))); 
			(ls(2,:).*sin(ls(1,:)))];
	lp2(:,:)=[ 	(ls(3,:).*cos(ls(1,:)+pi/2 )); 
			(ls(3,:).*sin(ls(1,:)+pi/2 ))];
case 3
	% We assume order (theta, phi, r1, r2, alpha2)
	% We now move a point 
	lp1(:,:)=[ 	(ls(3,:).*sin(ls(2,:)).*cos(ls(1,:))); 
			(ls(3,:).*sin(ls(2,:)).*sin(ls(1,:))); 
			(ls(3,:).*cos(ls(2,:)))];
	lp2(:,:)=[ 	(ls(4,:).*sin(ls(2,:)).*cos(ls(1,:)+pi/2)); 
			(ls(4,:).*sin(ls(2,:)).*sin(ls(1,:)+pi/2)); 
			(ls(4,:).*cos(ls(2,:)))];
	%lp3(:,:)=[ (ls(3,:).*sin(ls(2,:)).*cos(ls(1,:))); (ls(3,:).*sin(ls(2,:)).*sin(ls(1,:))); (ls(3,:).*cos(ls(2,:)))];
endswitch

% Plot some stuff for debugging purposes, only use if you have a few clusters
if (plot_avg_cov)
	printf("The averages are:\n")
	mu 
	printf("And the covariance matrices:\n")
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
	case "normal"
		% points on line segment are generated normally
		sn=randn(1,2);
	case "uniform"
		% points on line segment are generated uniformly
		sn=2*rand(1,2)-1;
	endswitch

	if (extremes_only) 
		% We set one of the coordinates to 1 or -1 with a certain probability
		pl1=(rand(1)>0.5)+1;
		pl2=2*(rand(1)>0.5)-1;
		sn(pl1)=pl2;
	endif

	% We use a multiplicative structure to pick a point on the line segment
	% The end of the line (lp) compared to the mean is taken negative as well as positive to form a line segment
	P(i,:)=mu(:,c) + lp1(:,c)*sn(1) + lp2(:,c)*sn(2) + chol(covar(:,:,c))' * randn(dim,1);
end

P=P';

switch (dim)
case 2
	% In the 2D case we plot a bit more, e.g. the "true" lines that generate the points as well as the mean values.
	plot(P(1,:),P(2,:),'.')
	return;
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



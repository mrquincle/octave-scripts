#!/usr/bin/octave -qf

% Not a function file
1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to generate N objects with a specified number of outliers and Gaussian noise
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
% Clear all previous sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear -x input_file output_file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Command line options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arg_list = argv ();

if (nargin > 1) 
	input_file = arg_list{1};
	output_file = arg_list{2};
else
	fprintf('Usage: dpmobjectrnd input output [skip_graphs]\n');
end

if (~exist("input_file"))
	fprintf('Error: Variable input_file not set\n');
	return
end

if (~exist("output_file"))
	fprintf('Error: Variable output_file not set\n');
	return
end

if (~exist(input_file))
	fprintf('Error: Input file does not exist\n');
	return
end

skip_graphs = false;

if (nargin == 3) 
	if (arg_list{3} == 'skip_graphs') 
		skip_graphs = true;
	end
end

graph_file=[output_file '.png'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default configuration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure Matlab-compatibility by casting broadcast warnings into errors
warning ('error', 'Octave:broadcast');

% add function to check if we run matlab or octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% include path for the Chinese Restaurant Process algorithm and the Normal-Inverse Wishart Random generator
addpath(genpath([pwd '/util']), 1);
addpath(genpath([pwd '/shapes']), 1);

if (~isOctave)
	% rehash, or changes to the configuration from config file will not be applied
	rehash
end

if (!skip_graphs)
set(gca(), "defaultlinelinewidth", 4.0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pkg load statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load configuration from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input_file=[config_dir '/' fname '_config.m'];
%output_file=[data_dir '/' fname '.pnts.data'];

if (~isOctave)
	% Matlab requires explicit clearing or else config file is not updated after a change
	clear(input_file);
end

fprintf('Load configuration file %s\n', input_file);
run(input_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print some configuration variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Object dimensions: %i\n', object_dim);
fprintf('Number of samples: %i\n', n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We run a Chinese Restaurant Process (CRP) for n items resulting a Dirichlet Process (DP) distribution. This is a 
% function that calculates the probabilities for all items in a batch. So, we run it only once. The resulting partition 
% is fixed.
[partitions clustersize clusters partition_bins] = crprnd(alpha,n);

% This means there is a fixed number of clusters rolling out of the CRP process above.
S=clusters;
fprintf('Number of clusters: %i\n', S);

% Initialize the parameters to the right dimensions
%  mu: mean
%  covar: covariance
%  lp: line endpoint
%  lp_covar: not used
%  P: points on the calculated line segment
mu=zeros(dim,S);
covar=zeros(dim,dim,S);

lss=zeros(param1.dim,S);

lp1=zeros(dim,S);
lp1_covar=zeros(dim,dim,S);
lp2=zeros(dim,S);
lp2_covar=zeros(dim,dim,S);
lp3=zeros(dim,S);
P=zeros(n,dim);
Y=zeros(n,dim);

% We add an additional CRP to get the same angle for line segments with a probability larger than 0, indicating 
% orrelations coming from higher-level objects, such as squares and alike.
[super_partitions super_clustersize super_clusters super_partition_bins] = crprnd(param1.alpha,S);
fprintf('Number of super clusters: %i\n', super_clusters);

% Normally distributed
s_mu=zeros(param1.dim,super_clusters);
s_covar=zeros(param1.dim,param1.dim,super_clusters);

% Pick the properties of the objects
for i=1:super_clusters
	[s_mu(:,i) s_covar(:,:,i)]=normalinvwishrnd(hyper1);
end

% We sample for each component from two different latent spaces. One corresponding to the center of a line segment, the
% other corresponding to its angle and length.

% For each cluster we get a new 'mean' vector and 'covariance' matrix for the center of the line segment.
% For each cluster we also get a new 'mean' vector and 'covariance' matrix for each (angle, segment size)
for i=1:S
	% We get a mean and covariance matrix at once through the normal-inverse Wishart random generator
	[mu(:,i) covar(:,:,i)]=normalinvwishrnd(hyper);

	% We get a mean and covariance matrix for the endpoint, but will throw away the covariance matrix later
	% This needs attention! The proper prior for an angle is different!
	%[lp(:,i) lp_covar(:,:,i)]=normalinvwishrnd(hyper);

	% Get the super-cluster
	super_cluster = super_partitions(i);

	% Get the angle and two lengths for 2D, get three angles and two lengths for 3D
	lss(:,i)=s_mu(:,super_cluster);

	% create rectangles, not parallelograms
	if (test_same_angle)
		switch(dim)
		case 2
			lss(1,:) = same_angle;
		case 3
			lss(1,:) = same_angle;
			lss(2,:) = same_angle;
		end
	end

	if (test_same_length) 
		switch(dim)
		case 2
			lss(2,:) = same_length;
			lss(3,:) = same_length;
		case 3
			lss(3,:) = same_length;
			lss(4,:) = same_length;
			lss(5,:) = same_length;
		end
	end
end

% For each pair of (angle,distance) calculate end point
switch(dim) 
case 2
	% We assume order (angle, length1, length2), we make a 'square' from it, actually a parallelogram by picking two
	% other points. Then all points will be between the middle \mu and points \mu+-lp1 and \mu+-lp2.
	lp1(:,:)=[ 	(lss(2,:).*cos(lss(1,:))); 
			(lss(2,:).*sin(lss(1,:)))];
	lp2(:,:)=[ 	(lss(3,:).*cos(lss(1,:)+pi/2 )); 
			(lss(3,:).*sin(lss(1,:)+pi/2 ))];
case 3
	switch (object_dim) 
	case 3
		% We assume order (theta, phi, r1, r2, r3, alpha2, alpha2)
		% We now move a point 
		lp1(:,:)=[ 	(lss(3,:).*sin(lss(2,:)).*cos(lss(1,:))); 
				(lss(3,:).*sin(lss(2,:)).*sin(lss(1,:))); 
				(lss(3,:).*cos(lss(2,:)))];
		lp2(:,:)=[ 	(lss(4,:).*sin(lss(2,:)).*cos(lss(1,:)+pi/2)); 
				(lss(4,:).*sin(lss(2,:)).*sin(lss(1,:)+pi/2)); 
				(lss(4,:).*cos(lss(2,:)))];
		lp3(:,:)=[ 	(lss(5,:).*sin(lss(2,:)+pi/2).*cos(lss(1,:))); 
				(lss(5,:).*sin(lss(2,:)+pi/2).*sin(lss(1,:))); 
				(lss(5,:).*cos(lss(2,:)+pi/2))];

	case 2
		% We assume order (theta, phi, r1, r2, alpha2)
		% We now move a point 
		lp1(:,:)=[ 	(lss(3,:).*sin(lss(2,:)).*cos(lss(1,:))); 
				(lss(3,:).*sin(lss(2,:)).*sin(lss(1,:))); 
				(lss(3,:).*cos(lss(2,:)))];
		lp2(:,:)=[ 	(lss(4,:).*sin(lss(2,:)).*cos(lss(1,:)+pi/2)); 
				(lss(4,:).*sin(lss(2,:)).*sin(lss(1,:)+pi/2)); 
				(lss(4,:).*cos(lss(2,:)))];
	end
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
		sn=randn(1,object_dim);
	case 'uniform'
		% points on line segment are generated uniformly
		% this doesn't really work nicely if the rectangles have lengths that differ in the x and y directions
		sn=2*rand(1,object_dim)-1;
	end

	if (extremes_only) 
		% We set all the coordinates to 1 or -1 with a certain probability except for one
		%pl1=(rand(1)>0.5)+1;

		% create select mask for which one will not be altered, say (0, 1, 0)
		plA=zeros(1,object_dim);
		plA(randi(object_dim))=1;
		% create select mask for the +1/-1 items, (1, 0, 1)
		plB=1-plA;
		% calculate +1/-1 items, so we get mask (-1, 0, +1)
		plB=plB.*((rand(1,object_dim)>0.5)*2-1);

		% apply both masks
		sn=sn.*plA+plB;
		%plA=unidrnd(dim);
		%plB=2*(rand(1)>0.5)-1;
		%sn(plA)=plB;
		%if (enable_cubes)
	end

	% We use a multiplicative structure to pick a point on the line segment
	% The end of the line (lp) compared to the mean is taken negative as well as positive to form a line segment
	switch (object_dim)
	case 2
		P(i,:)=mu(:,c) + lp1(:,c)*sn(1) + lp2(:,c)*sn(2) + chol(covar(:,:,c))' * randn(dim,1);
	case 3
		% Just forget all of the above, just sample cubes...
		e = sample_cube(lss(3,c),lss(4,c),lss(5,c)) - [lss(3,c) lss(4,c) lss(5,c)]/2;
		P(i,:) = sample_line(e)' + mu(:,c) + chol(covar(:,:,c))' * randn(dim,1);
	%	P(i,:)=mu(:,c) + lp1(:,c)*sn(1) + lp2(:,c)*sn(2) + lp3(:,c)*sn(3) + chol(covar(:,:,c))' * randn(dim,1);
	end
	Y(i,:)=c;
end

% Transpose matrix for easy plotting
P=P';
Y=Y';

% Scale points
scale=100;
P=P./scale;
mu=mu./scale;

switch (dim)
case 2
	dlmwrite(output_file, [P(1,:)' P(2,:)', Y(1,:)'], 'delimiter', '\t', 'precision', 10);
case 3
	dlmwrite(output_file, [P(1,:)' P(2,:)' P(3,:)', Y(1,:)'], 'delimiter', '\t', 'precision', 10);
end

if (!skip_graphs)
	switch (dim)
	case 2
		% In the 2D case we plot a bit more, e.g. the 'true' lines that generate the points as well as the mean values.
		plot(P(1,:),P(2,:),'.', 'MarkerSize', 12)
		dlmwrite(output_file, [P(1,:)' P(2,:)', Y(1,:)'], 'delimiter', '\t', 'precision', 10);
		hold on
		% Print mean values as a red circle
		plot(mu(1,:),mu(2,:),'or');
	case 3
		plot3(P(1,:),P(2,:),P(3,:),'.')
		hold on
		% Print mean values as a red circle
		plot3(mu(1,:),mu(2,:),mu(3,:),'or');
		dlmwrite(output_file, [P(1,:)' P(2,:)' P(3,:)', Y(1,:)'], 'delimiter', '\t', 'precision', 10);
	end

	axis equal

	print(graph_file, "-dpng")

end

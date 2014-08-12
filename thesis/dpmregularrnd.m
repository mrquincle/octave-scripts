#!/usr/bin/octave -qf
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
% Goal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This specific file has as goal to come up with some method to introduce regular spacing as often seen in man-made
% environments. For example the windows in a flat are regularly spaced. This structure exists in the "head of the 
% designer" and reflect engineering principles of building robustly and affordably.

% The way we incorporate this concept is by introducing yet another layer of "spacing objects", which can be seen as
% the facade of a building for example, but calling it like that is just for our convenience. The inference method
% will just have a way to collect spatially regular items. 

% We first choose a "type" and then we choose a "spacing object". This means that different types can still intersect
% each other. And even intersection of spacing objects of the same type is possible.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load configuration from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname="dpmregular1";
config_dir="config";

data_dir="data";
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate items from a CRP to model a statistic process on hierarchy level 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We run a Chinese Restaurant Process (CRP) for n items resulting a Dirichlet Process (DP) distribution. This is a 
% function that calculates the probabilities for all items in a batch. So, we run it only once. The resulting partition 
% is fixed.
[h0_p h0_cs h0_c h0_pb] = crprnd(hyper0.alpha,n);

% This means there is a fixed number of clusters rolling out of the CRP process above.
printf("Number of clusters on level 0: %i\n", h0_c);

% Initialize the parameters to the right dimensions
%  mu: mean
%  covar: covariance
%  lp: line endpoint
%  lp_covar: not used
%  P: points on the calculated line segment
mu=zeros(hyper0.dim,h0_c);
covar=zeros(hyper0.dim,hyper0.dim,h0_c);

ls=zeros(hyper1.dim,h0_c);
lg=zeros(hyper2.dim,h0_c);

lp1=zeros(hyper0.dim,h0_c);
lp1_covar=zeros(hyper0.dim,hyper0.dim,h0_c);
lp2=zeros(hyper0.dim,h0_c);
lp2_covar=zeros(hyper0.dim,hyper0.dim,h0_c);
lp3=zeros(hyper0.dim,h0_c);
P=zeros(n,hyper0.dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate items from a CRP to model a statistic process on hierarchy level 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We add an additional CRP to get the same angle for line segments with a probability larger than 0, indicating 
% correlations coming from higher-level objects, such as squares and alike.
[h1_p h1_cs h1_c h1_pb] = crprnd(hyper1.alpha, h0_c);
printf("Number of clusters on level 1: %i\n", h1_c);

% Normally distributed
h1_mu=zeros(hyper1.dim, h1_c);
h1_covar=zeros(hyper1.dim, hyper1.dim, h1_c);

% Pick the properties of the objects
for i=1:h1_c
	[h1_mu(:,i) h1_covar(:,:,i)]=normalinvwishrnd(hyper1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate items from a CRP to model a statistic process on hierarchy level 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h2_p h2_cs h2_c h2_pb] = crprnd(hyper2.alpha, h1_c);
printf("Number of clusters on level 2: %i\n", h2_c);

% Normally distributed
h2_mu=zeros(hyper2.dim, h2_c);
h2_covar=zeros(hyper2.dim, hyper2.dim, h2_c);

% Pick the properties of the objects on this level
for i=1:h2_c
	[h2_mu(:,i) h2_covar(:,:,i)]=normalinvwishrnd(hyper2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We sample for each component from two different latent spaces. One corresponding to the center of an object, the
% other corresponding to its properties. 

% In the case of a line segment, this is the orientation of the line and the size of the segment. In the case of a
% parallelogram the vector "ls" should be larger and contains additional parameters to specify such a more complex
% object.

for i=1:h0_c
	% We get a mean and covariance matrix at once through the normal-inverse Wishart random generator
	[mu(:,i) covar(:,:,i)]=normalinvwishrnd(hyper0);

	% Get the cluster id
	h1_ci = h1_p(i);
	% Get the necessary parameter values for our object (same for all objects in this cluster)
	ls(:,i)=h1_mu(:,h1_ci);

	% Get the grid parameters
	h2_ci = h2_p(h1_ci);
	lg(:,i)=h2_mu(:,h2_ci);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate properties of the object using the generated random variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate important points of the object to be drawn. For a line segment this is a single point, which is subsequentl
% used together with the mean \mu to draw a line segment. The other end of the segment does not need to be calculated 
% because the points will be drawn from \mu-lp1 to \mu+lp1 (so the point is "mirrored" over the mean point). Hence, the
% point is not an absolute position, but relative. This makes it easier to create objects from the same size.
% The other objects, parallelograms, cubes, etc. are created similarly.
switch(hyper0.dim) 
case 2
	switch (object_dim) 
	case 1
		% Assumed order in "ls" vector: (angle, distance)
		lp1(:,:)=[ 	(ls(2,:).*cos(ls(1,:))); 
				(ls(2,:).*sin(ls(1,:)))];
	case 2
		% We assume order (angle, length1, length2), we make a "square" from it, actually a parallelogram by picking two
		% other points. Then all points will be between the middle \mu and points \mu+-lp1 and \mu+-lp2.
		lp1(:,:)=[ 	(ls(2,:).*cos(ls(1,:))); 
				(ls(2,:).*sin(ls(1,:)))];
		lp2(:,:)=[ 	(ls(3,:).*cos(ls(1,:)+pi/2 )); 
				(ls(3,:).*sin(ls(1,:)+pi/2 ))];
	endswitch
case 3
	switch (object_dim) 
	case 2
		% We assume order (theta, phi, r1, r2, alpha2)
		lp1(:,:)=[ 	(ls(3,:).*sin(ls(2,:)).*cos(ls(1,:))); 
				(ls(3,:).*sin(ls(2,:)).*sin(ls(1,:))); 
				(ls(3,:).*cos(ls(2,:)))];
		lp2(:,:)=[ 	(ls(4,:).*sin(ls(2,:)).*cos(ls(1,:)+pi/2)); 
				(ls(4,:).*sin(ls(2,:)).*sin(ls(1,:)+pi/2)); 
				(ls(4,:).*cos(ls(2,:)))];
	case 3
		% We assume order (theta, phi, r1, r2, r3, alpha2, alpha2)
		lp1(:,:)=[ 	(ls(3,:).*sin(ls(2,:)).*cos(ls(1,:))); 
				(ls(3,:).*sin(ls(2,:)).*sin(ls(1,:))); 
				(ls(3,:).*cos(ls(2,:)))];
		lp2(:,:)=[ 	(ls(4,:).*sin(ls(2,:)).*cos(ls(1,:)+pi/2)); 
				(ls(4,:).*sin(ls(2,:)).*sin(ls(1,:)+pi/2)); 
				(ls(4,:).*cos(ls(2,:)))];
		lp3(:,:)=[ 	(ls(5,:).*sin(ls(2,:)+pi/2).*cos(ls(1,:))); 
				(ls(5,:).*sin(ls(2,:)+pi/2).*sin(ls(1,:))); 
				(ls(5,:).*cos(ls(2,:)+pi/2))];
	endswitch
endswitch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate properties of grid using the random variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor_spacing=20;
gdx=lg(1,:)*factor_spacing;
gdy=lg(2,:)*factor_spacing;
gnx=floor(abs(lg(3,:))+1);
gny=floor(abs(lg(4,:))+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print some interesting information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printf("There are %i objects of %i types in %i grids\n", h0_c, h1_c, h2_c);
for i=1:h0_c
	h1_ci = h1_p(i);
	h2_ci = h2_p(h1_ci);
	printf("Object %i of type %i \n", i, h1_ci);
	printf("Location: ");
	disp(mu(:,i)');
	printf("Object parameters: ");
	disp(h1_mu(:,h1_ci)');
	switch (object_dim) 
	case 1
		obj_size=abs(ls(2,i));
		printf("Object size: %d\n", obj_size);
	case 2
		obj_size1=abs(ls(2,i));
		obj_size2=abs(ls(3,i));
		printf("Object size: %fx%f\n", obj_size1, obj_size2);
	endswitch
	printf("Grid layout: ");
	disp(h2_mu(:,h2_ci)');
	printf("This means a grid of %ix%i with spacing %f and %f\n", gnx(i), gny(i), gdx(i), gdy(i));
endfor

% I hate to press q or enter all the time
fflush(stdout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print some other stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
	c=h0_p(i);

	% Line distribution
	switch (noise_distribution) 
	case "normal"
		% points on line segment are generated normally
		sn=randn(1,object_dim);
	case "uniform"
		% points on line segment are generated uniformly
		sn=2*rand(1,object_dim)-1;
	endswitch

	if (extremes_only) 
		% We set all the coordinates to 1 or -1 with a certain probability except for one

		% create select mask for which one will not be altered, say (0, 1, 0)
		plA=zeros(1,object_dim);
		plA(unidrnd(object_dim))=1;
		% create select mask for the +1/-1 items, (1, 0, 1)
		plB=1-plA;
		% calculate +1/-1 items, so we get mask (-1, 0, +1)
		plB=plB.*((rand(1,object_dim)>0.5)*2-1);

		% apply both masks
		sn=sn.*plA+plB;
	endif

	% Transform grid parameters in shift operator for points
	nxi=unidrnd(gnx(c));
	nyi=unidrnd(gny(c));
	gp(:,c) = [gdx(c) * nxi; gdy(c) * nyi ];


	% We use a multiplicative structure to pick a point on the line segment
	% The end of the line (lp) compared to the mean is taken negative as well as positive to form a line segment
	switch (object_dim)
	case 1
		P(i,:)=mu(:,c) + lp1(:,c)*sn + chol(covar(:,:,c))' * randn(hyper0.dim,1);
	case 2
		P(i,:)=mu(:,c) + gp(:,c) + lp1(:,c)*sn(1) + lp2(:,c)*sn(2) + chol(covar(:,:,c))' * randn(hyper0.dim,1);
	case 3
		P(i,:)=mu(:,c) + lp1(:,c)*sn(1) + lp2(:,c)*sn(2) + lp3(:,c)*sn(3) + chol(covar(:,:,c))' * randn(hyper0.dim,1);
	endswitch
end

P=P';

fh=figure(1);

show_arrows=false;

for i=1:1
switch (hyper0.dim)
case 2
	disp("Show graph");
	% In the 2D case we plot a bit more, e.g. the "true" lines that generate the points as well as the mean values.
	plot(P(1,:),P(2,:),'.')
	if (!show_arrows) 
		disp("Break out of loop. Do not show arrows and print to file");
		break;
	endif
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
	for i=1:h0_c
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
endfor

disp("Wait till plot has been closed");
uiwait(fh);
disp("Thanks for everything!");
fflush(stdout);

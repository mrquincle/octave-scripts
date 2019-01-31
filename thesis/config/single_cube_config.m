%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration file for dpmobjectrnd.m
%
% Author: Anne van Rossum
% Date: Jul, 2014
% Copyrights: Almende B.V., Rotterdam, The Netherlands
% License: LGPLv3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_avg_cov = false;

% A 3D world
dim=3;
% A line is a 2D object, set to 2 in that case!
object_dim=3;
% the number of items is fixed, it is not possible to get another item in an incremental fashion using this generator
n=2048;
% With a higher alpha there will be more cubes formed: the chance to form a new cluster will increase
% For now set alpha almost zero, so we have only a single cube.
alpha=0.000001;

% distribution of objects over the entire space

% mean for each cluster, D=2
hyper.mu = zeros(dim,1);
% used to multiply the covariance matrix with
hyper.kappa = 0.0001;
% nu should be larger than D-1
hyper.nu = dim+2;
% covariance matrix
% the spread of the objects over the entire space
hyper.lambda = 8*eye(dim);

% A larger value for param1.alpha means more different types of lines
% The dimension is now (param1.alpha, length1, length2) for 2D, and (param1.alpha, alpha2, length1, length2) for 3D
switch (dim)
case 2
	param1.dim=dim+1;
case 3
	param1.dim=dim+2;
end

switch (object_dim)
case 3
	if (dim == 2) 
		error ('It does not make sense to try to generate 3D objects in a 2D world')
	end
	param1.dim = param1.dim + 1;
end

% distribution of points over an object

param1.alpha=1;

hyper1.mu = zeros(param1.dim,1);
% variance of the points on the lines / squares
hyper1.kappa = 1;
hyper1.nu = param1.dim+1;
% size of the object
hyper1.lambda = eye(param1.dim) * 1;

% choose noise distribution to be 'normal' or 'uniform'
%noise_distribution='normal';
noise_distribution='uniform';

test_same_angle=true;
same_angle=0;

test_same_length=true;
same_length=100;

extremes_only=true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration file for dpmobject.m
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

dim=3;
% With a higher alpha there will be more lines formed, the chance to form a new cluster will increase
alpha=2;

% mean for each cluster, D=2
hyper.mu = zeros(dim,1);
% used to multiply the covariance matrix with
hyper.kappa = 0.001;
% nu should be larger than D-1
hyper.nu = dim+2;
% covariance matrix
% the spread of the objects over the entire space
hyper.lambda = 8*eye(dim);

% the number of items is fixed, it is not possible to get another item in an incremental fashion using this generator
n=2000;

% A larger value for param1.alpha means more different types of lines
% The dimension is now (param1.alpha, length1, length2) for 2D, and (param1.alpha, alpha2, length1, length2) for 3D
switch (dim)
case 2
	param1.dim=dim+1;
case 3
	param1.dim=dim+2;
end

% A line is a 2D object, set to 2 in that case!
object_dim=3;
switch (object_dim)
case 3
	if (dim == 2) 
		error ('It does not make sense to try to generate 3D objects in a 2D world')
	end
	param1.dim = param1.dim + 1;
end

param1.alpha=1;

hyper1.mu = zeros(param1.dim,1);
hyper1.kappa = 0.001;
hyper1.nu = param1.dim+1;
hyper1.lambda = 1*eye(param1.dim);

% choose noise distribution to be 'normal' or 'uniform'
%noise_distribution='normal';
noise_distribution='uniform';

test_same_angle=false;
same_angle=pi/3;

test_same_length=false;
same_length=4;

extremes_only=true;

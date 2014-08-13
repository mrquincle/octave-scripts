%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration file for dpmline.m
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

dim=2;
alpha=4;

% mean for each cluster, D=2
hyper.mu = zeros(dim,1);
% used to multiply the covariance matrix with
hyper.kappa = 0.0001;
% nu should be larger than D-1
hyper.nu = 4;
% covariance matrix
hyper.lambda = 1*eye(dim);


% the number of items is fixed, it is not possible to get another item in an incremental fashion using this generator
n=200;

% choose noise distribution to be 'normal' or 'uniform'
noise_distribution='normal';

test_same_angle=false;
same_angle=pi/3;

test_same_length=false;
same_length=4;

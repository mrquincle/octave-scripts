%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration file for onedim.m
%
% Author: Anne van Rossum
% Date: Feb, 2015
% Copyrights: Almende B.V., Rotterdam, The Netherlands
% License: LGPLv3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_avg_cov = false;

dim=1;
% choose alpha larger for more clusters (alpha > 0)
alpha=0.7;

% mean for each cluster, D=2
hyperA.mu = zeros(dim,1);
% used to multiply the covariance matrix with
hyperA.kappa = 0.005;
% nu should be larger than D-1
hyperA.nu = 4;
% covariance matrix
hyperA.lambda = 1*eye(dim);

% mean for each cluster, D=2
hyperB.mu = zeros(dim,1);
% used to multiply the covariance matrix with
hyperB.kappa = 1;
% nu should be larger than D-1
hyperB.nu = 4;
% covariance matrix
hyperB.lambda = 1*eye(dim);

% the number of items is fixed, it is not possible to get another item in an incremental fashion using this generator
n=200;

% choose noise distribution to be 'normal' or 'uniform'
noise_distribution='unform';
% should actually also come from a prior
sn_scale = 0.01;
sn_scale = 1;

test_same_angle=false;
same_angle=pi/3;

test_same_length=false;
same_length=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to generate uniform distributions in one dimenion
%
% Author: Anne van Rossum
% Date: Feb, 2015
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

fname='onedim';
config_dir='config';

data_dir='data';
mkdir(data_dir);

input_file=[config_dir '/' fname '_config.m']
output_file=[data_dir '/' fname '.pnts.data'];

run(input_file);

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

% generate partitions
[partitions clustersize clusters partition_bins] = crprnd(alpha,n);

% Initialize the parameters to the right dimensions
%  S: # of clusters
%  mu: mean
%  covar: covariance
%  lp: line endpoint
%  lp_covar: not used
%  P: points on the calculated line segment
%  Y: assignments
S=clusters;
mu=zeros(dim,S);
covar=zeros(dim,dim,S);
lp=zeros(dim,S);
lp_covar=zeros(dim,dim,S);
P=zeros(n,dim);
Y=zeros(n,dim);

% for each cluster (line segment) generate two end points
for i=1:S
	[mu(:,i) covar(:,:,i)]=normalinvwishrnd(hyperA);

	[lp(:,i) lp_covar(:,:,i)]=normalinvwishrnd(hyperB);
end

% generate points in between the end points
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
	sn=sn*sn_scale;

        % We use a multiplicative structure to pick a point on the line segment
        % The end of the line (lp) compared to the mean is taken negative as well as positive to form a line segment
        P(i,:)=mu(:,c) + lp(:,c)*sn + chol(covar(:,:,c))' * randn(dim,1);
        Y(i,:)=c;
end

% plot
for i = 0:S
	yp=[P,Y];
	yp(find(yp(:,2)!=i),:)=[];
	clr=hsv2rgb([i/S,1,1]);
	plot(yp(:,1),0,'color',clr,'.');
	hold on
end

dlmwrite(output_file, [P,Y], '\t', "precision", 10);

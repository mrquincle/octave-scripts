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

fname="dpmgauss1";
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

% run Dirichlet process for n items
% we now have a DP distribution
% note, that we have a fixed partition
[partitions clustersize clusters partition_bins] = crprnd(alpha,n);

S=clusters;
printf("Number of clusters: %i\n", clusters);

mu=repmat(zeros(1,2),S,1)';
covar=repmat(zeros(2,2),S,1)';

% loop to get S items with a mean and covariance matrix
for i=1:S
	%[mu(i) samples(i)]=
	% get combined mean and covariance matrix
	[mu(:,i) covar(:,i*2-1:i*2)]=normalinvwishrnd(hyper);
end

if (plot_avg_cov)
	printf("The averages are:\n")
	mu 
	printf("And the covariance matrices:\n")
	covar
end

P=zeros(n,2);
for i=1:n
	% Sample assignment or cluster index from CRP
	c=partitions(i);
	% Sample point from mean 
	P(i,:)=mu(:,c) + chol(covar(:,c*2-1:c*2))' * randn(length(mu(:,c)),1);
end

P=P';
plot(P(1,:),P(2,:),'.')

hold on

% print mean values
plot(mu(1,:),mu(2,:),'or');

dlmwrite(output_file, [P(1,:)' P(2,:)'], '\t', "precision", 10);

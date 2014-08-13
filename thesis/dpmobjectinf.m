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

fname='dpmobject1';
config_dir='config';

data_dir='data';

input_file=[config_dir '/' fname '_config.m'];
data_file=[data_dir '/' fname '.pnts.data'];

run(input_file);

% Load data in P
P=load(data_file)';

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

% hyper and alpha come from config file
hyper1
alpha=0.1
niter=100;
doPlot=2;
type_algo='collapsedCRP';
[c_st, c_est, similarity] = gibbsDPM(P, hyper1, alpha, niter, type_algo, doPlot);

[~, ind] = sort(c_est);
figure
imagesc(double(similarity(ind, ind)));

figure
ind2 = unique(c_est);
for j=1:length(ind2)
	subplot(5,5,j)
	plot(P(:, c_est==ind2(j)), 'k')
	ylim([-5, 5])
end


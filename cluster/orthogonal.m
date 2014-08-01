%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of samples
n=10;

show_plots=true;
config_use_noise=1;
test_introduce_outliers=1;
test_outliers_low_weight=1; % if outliers have the same weight, they screw up everything!
compare_deming=1;

% ensure Matlab-compatibility by making broadcast warnings, errors
warning ("error", "Octave:broadcast");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create sample points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate (x,y) random pairs with white noise in both x and y direction

% line coefficients
g0=2;
g1=1;
printf("Test with line: y=%f+%f*x\n", g0, g1);

% create n points on this line
r=rand(1,n)/n;
if (!config_use_noise) 
	r=zeros(1,n);
endif
x=[1:n]/n + r;
r=rand(1,n)/n;
if (!config_use_noise) 
	r=zeros(1,n);
endif
y=g0+g1*x + r;
pnts=[x; y]';
x_axis_len=n/n;

% introduce two outliers
if (test_introduce_outliers)
%	pnts(1,:)=[.5, 6];
	pnts(1,:)=[1.5, 3];
	pnts(2,:)=[4,15];
endif

x=pnts(:,1);
y=pnts(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare with function from the web
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (compare_deming)
	dem=deming(x,y);
	printf("Deming according to function from the web: y=%f+%f*x\n", dem(1), dem(2));
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weight the points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for now weight everything just as much, except point 1 and 2
weights=ones(n,1)*1;
if (test_outliers_low_weight)
	weights(1)=0.1;
	weights(2)=0.1;
endif
total_weight=sum(weights);

% calculate weighted points
wpnts = pnts.*repmat(weights,1,2);

% calculate average of weighted points
avg=sum(wpnts) / total_weight;

avg_x=avg(1);
avg_y=avg(2);
printf("Sample averages: M[x]=%f, M[y]=%f\n", avg_x, avg_y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move all points with the sample average, so they are centered around the origin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate x_i - x_avg and y_i - y_avg
diff=bsxfun(@minus, pnts, avg);

% our centered points are our original unweighted points, but shifted with the weighted average
cpnts=diff;
%cpnts=pnts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recalculate weights and adjust position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjusted points, will be centered and weighted
apnts = cpnts.*repmat(weights,1,2);

% https://en.wikipedia.org/wiki/Moment_(mathematics)
% raw sampling is divide / n
% second central moment: 1 / (n-1)
moment_factor=(n-1);

% adjusted sample variance 
var=sum(apnts.^2)/moment_factor; 
sxx=var(1); syy=var(2);
sxy=sum(apnts(:,1).*apnts(:,2))/moment_factor;
printf("Sample variances: S[x,x]=%f, S[y,y]=%f, S[x,y]=%f\n", sxx, syy, sxy);

sgn=1;
lambda=1;
b1=(syy-lambda*sxx+sgn*sqrt( (syy-lambda*sxx).^2 + 4*lambda*sxy.^2)) / (2*sxy);
b0=avg_y - b1*avg_x;

aavg=sum(apnts);
aavg_x=aavg(1);
aavg_y=aavg(2);
printf("Shifted sample average should be (0,0): M[x]=%f, M[y]=%f\n", aavg_x, aavg_y);

% in case you want to plot the line through the origin
c1=(syy-sxx+sqrt( (syy-sxx)^2 + 4*sxy^2)) / (2*sxy);
c0=aavg_y - c1*aavg_x;

% show coefficients of line to user
printf("Calculated line: y=%f+%f*x\n", b0, b1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (show_plots) 
	subplot(2,2,1)
	text(.75, 1.25, 'Orthogonal regression under outliers')

	% plot the points themselves
	plot(pnts(:,1), pnts(:,2), '*5');
	%plot(pnts(:,1), pnts(:,2), 'ro', "markersize",10, "linewidth",3);
	hold on;
	
	% plot the average
	navg=sum(pnts)/size(pnts,1);
	plot(navg(1), navg(2), 'o1');
	
	% draw the estimated line 
	line=[0, b0; x_axis_len, b0+b1*x_axis_len]';
	plot(line(1,:), line(2,:), '-');
	title('Original points');
	axis('square');
	axis('equal');
	axis('auto');
	xlabel('Variable x with Gaussian noise');
	ylabel('Variable y with Gaussian noise');
	legend('(x,y) pairs', 'average', 'regression line');
	legend location east
	
	hold off;

	subplot(2,2,2)

	% plot the weighted points themselves
	plot(wpnts(:,1), wpnts(:,2), '*5');
	hold on;
	
	% plot the average
	wavg=sum(wpnts) / total_weight;
	plot(wavg(1), wavg(2), 'o1');
	
	% draw the estimated line 
	line=[0, b0; x_axis_len, b0+b1*x_axis_len]';
	plot(line(1,:), line(2,:), '-');
	title('Points with outliers with lower weights');
	axis('square');
	axis('equal');
	axis('auto');
	xlabel('Variable x with Gaussian noise');
	ylabel('Variable y with Gaussian noise');
	legend('(x,y) pairs', 'average', 'regression line');
	legend location east
	
	hold off;

	subplot(2,2,3)
	plot(cpnts(:,1), cpnts(:,2), '*5');

	hold on;
	line=[-x_axis_len, -b1*x_axis_len; x_axis_len, b1*x_axis_len]';
	plot(line(1,:), line(2,:), '-');
	
	title('Points shifted with the weighted average')

	axis('square');
	axis('equal');
	xlabel('Variable x with Gaussian noise');
	ylabel('Variable y with Gaussian noise');
	legend('(x,y) pairs', 'regression line');
	legend location east

	hold off
	
	subplot(2,2,4)
	plot(apnts(:,1), apnts(:,2), '*5');

	hold on;
	line=[-x_axis_len, -b1*x_axis_len; x_axis_len, b1*x_axis_len]';
	plot(line(1,:), line(2,:), '-');
	
	title('Adjusted points, shifted and weighted')
	axis('square');
	axis('equal');
	xlabel('Variable x with Gaussian noise');
	ylabel('Variable y with Gaussian noise');
	legend('(x,y) pairs', 'regression line');
	legend location southeast
	hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating difference with original line description 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note that ordinary least squares does calculate the distance along one dimension (e.g. y) and it is not an 
% errors-in-variables model, with a distance calculated orthogonal to the line
dd=((b0+x*b1)-y).^2';
dist=sum(dd);
printf("Distance calculated in ordinary least squares (over y-axis): %f\n", dist);

% in case of total squares, in |ax+by+c| notation a=b1, b=-1, c=b0
%dd=(abs(x*b1-y+b0) / sqrt(b1.^2 + 1)).^2
dd=(x*b1-y+b0).^2 / (b1.^2 + 1);
dist=sum(dd);
printf("Distance calculated in total least squares: %f\n", dist);

% and in case of original line it would be
dd=((g0+x*g1)-y).^2';
dist=sum(dd);
printf("Distance calculated if original line would have been known (without outliers), in ordinary least squares: %f\n", dist);
dd=(x*g1-y+g0).^2 / (g1.^2 + 1);
dist=sum(dd);
printf("Distance calculated if original line would have been known (without outliers), in total least squares: %f\n", dist);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configuration parameters, from adding Gaussian noise, to test options
config_use_noise=0;
test_horizontal_line=0;
test_vertical_line=0;
test_calculate_first_two_points=0;
test_calculate_angle_first_two_points=0;
print_parameters=1;

test_line_estimation=1;
test_introduce_outlier=1;
do_threshold_weights=1;
plot_points=0;
plot_distance_angle=1;
plot_weighted_points=1;
plot_shifted_weighted_points=1;

% number of samples
n=40;

% the level of noise can be reduced or increased
noise_factor=1;

% ensure Matlab-compatibility by casting broadcast warnings into errors
warning ("error", "Octave:broadcast");

% line coefficients, format: [scale, y intersect; scale, y intersect], etc.
G=[1 2; -1 -1; -10 4];
G=[1 2; -1 -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate points according to above configuration settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensionality is normally two
dim=size(G,2);

% L is the number of lines
L=size(G,1);

% generate x coordinates
R=rand(L,n)/n*noise_factor;
if (!config_use_noise) 
	R=zeros(L,n);
endif
X=repmat([1:n]/n, L, 1) + R;

% generate y coordinates
R=rand(L,n)/n*noise_factor;
if (!config_use_noise) 
	R=zeros(L,n);
endif
Y=repmat(G(:,1), 1, n) .*X + repmat(G(:,2), 1, n) + R;

% generate cluster index
K=ceil([1/n:1/n:L])';

% result of the points with the line index, this can be written to a file
D=[X'(:) Y'(:) K];

% get only the points, without the line index
F=[X'(:) Y'(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate points according to above configuration settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% testing: if we have two points that are exactly vertically or horizontally oriented from each other
if (test_vertical_line) 
	F(2,:)=F(1,:);
	% introduce a vertical line by positioning the second point exactly above the first point
	F(2,2)=F(2,2)+1;
endif
if (test_horizontal_line)
	F(3,:)=F(1,:);
	% introduce a horizontal line
	F(3,1)=F(3,1)+1;
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a matrix with difference between points, and use it calculate the angle of the line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate a 3D array with the difference between points with each other
Fdiff=bsxfun(@minus, F, permute(F, [3 2 1]));

% calculate the angle as atan2(y,x)
A=atan2(Fdiff(:,2,:),Fdiff(:,1,:));

% modulus pi, so everything is nicely between 0 and pi
A=mod(A+2*pi, pi);

% reshape 
A=reshape(A,n*L,n*L)';

if (test_calculate_first_two_points)
	% rotate with second angle, rotate first point
	xnew=F(1,1)*cos(A(1,2))+F(1,2)*sin(A(1,2));
	ynew=-F(1,1)*sin(A(1,2))+F(1,2)*cos(A(1,2));
	printf("Rotate point with angle of %f\n", A(1,2));
	P1=F(1,:)
	N1=[xnew ynew]

	% rotate with second angle, rotate second point
	xnew=F(2,1)*cos(A(1,2))+F(2,2)*sin(A(1,2));
	ynew=-F(2,1)*sin(A(1,2))+F(2,2)*cos(A(1,2));
	printf("Rotate point with angle of %f\n", A(1,2));
	P2=F(2,:)
	N2=[xnew ynew]

	% calculate the distance
	x1=P1(1,1)
	y1=P1(1,2)
	x2=P2(1,1)
	y2=P2(1,2)
	dd=(-x1*y2+x2*y1)/sqrt((x2-x1)^2 + (y2-y1)^2 )
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the calculate angle to rotate the line and read off the distance to the origin by the y-intersection value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rotate all points around the origin
s=sin(A);
c=cos(A);
x_times_sin=-repmat(F(:,1),1,L*n).*s;
y_times_cos=repmat(F(:,2),1,L*n).*c;

% calculate the distance to the origin
D=abs(x_times_sin+y_times_cos).*(1-eye(L*n));

if (test_calculate_angle_first_two_points)
	% manually calculate pnt2-pnt1 and pnt3-pnt1 for testing
	dpnt=[pnt2-pnt1; pnt3-pnt1];
	ang=mod(atan2(dpnt(:,2), dpnt(:,1))+2*pi, pi);
endif

printf("Print angle and distance as pairs to file\n")
DA=[squareform(D) squareform(A)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show the parameters
if (print_parameters)
	beta=G
	avg=[mean(X'); mean(Y')]
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dlmwrite('lines.data', D, '\t', "precision", 10);
dlmwrite('linepairs.data', DA, '\t', "precision", 10);

if (plot_points)
	subplot(1,2,1)
	text(.75, 1.25, 'How noise shows up in distance-angle space')
	plot(X',Y','*');
	title('Original points with noise')
	lab={1:L};
	for i = [1:L]
		lab{i}=sprintf('line y=%i*x+%i', G(i,1), G(i,2));
	endfor

	legend(lab)
	%legend boxoff
	%legend location east
	xlabel('some coordinate')
	ylabel('some other coordinate')
	axis('square')
	axis('equal')
	grid on

	subplot(1,2,2)
	plot(DA(:,1),DA(:,2),'*')
	title('In distance-angle space')
	xlabel('distance to origin')
	ylabel('angle with origin')
	axis('square')
	axis('equal')
	
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform line estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The line estimation test can be done if you define above 2 clusters, with half of the points assigned to one, the 
% other half to the other cluster. 
if (test_line_estimation) 

	% first calculate weights, say we are already close to the first cluster
	weights=[ones(1,n)*0.9, ones(1,n)*0.1 ]';

	% quadratic weights is the same as repeating the line estimation step
	%weights=weights.^2;

	% introduce an outlier by assigning the 14th item a large weight
	if (test_introduce_outlier)
		weights(14)=0.3;
	endif

	if (do_threshold_weights)
		% if weight below 0.2, set weight to 0:
		id_w_thr=(weights < 0.2);
		weights(id_w_thr)=0;
	endif

	% weigh points 
	apnts = F.*repmat(weights,1,2);

	% calculate weighted average
	avg=sum(apnts)./sum(weights)

	% we move all points with the weighted average, so that if correct, the line is now oriented around the origin
	diff=bsxfun(@minus, F, avg);

	% we multiply the shifted points with their weights, this expands points on the line, and shrinks the others
	dpnts=diff.*repmat(weights,1,2);

	if (plot_weighted_points)
		plot(apnts(:,1), apnts(:,2), '+')
		hold on;
	endif

	if (plot_shifted_weighted_points)
		plot(dpnts(:,1), dpnts(:,2), '*')
		hold off;
	endif

	diff=dpnts;

	% https://en.wikipedia.org/wiki/Moment_(mathematics)
	% raw sampling is divide / n, second central moment: 1 / (n-1), we need the latter
	moment_factor=(n-1);

	% adjusted sample variance 
	var=sum(diff.^2)/moment_factor; 
	sxx=var(1); syy=var(2);
	sxy=sum(diff(:,1).*diff(:,2))/moment_factor;
	printf("Sample variances: S[x,x]=%f, S[y,y]=%f, S[x,y]=%f\n", sxx, syy, sxy);

	% calculate the line parameters b0 and b1
	lambda=1;
	sgn=1;
	b1=(syy-lambda*sxx+sgn*sqrt( (syy-lambda*sxx).^2 + 4*lambda*sxy.^2)) / (2*sxy);
	b0=avg(1) - b1*avg(2);
	printf("Estimated line coordinate: y=%f*x+%f\n", b1, b0);

	% TODO: check, it seems refactoring introduced a negative sign in b0
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% It might be interesting to observe a histogram. Take all the points and calculate their distance to the line. Then
% create a histogram in which points are divided over distance intervals.
%hist(dpnts, 30);


% This is especially interesting with respect to "wrong" lines. How do the residuals show up, and more important, how 
% to move towards a better fit!?




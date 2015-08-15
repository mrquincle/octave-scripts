% two-sided Uniform likelihood
% we assume lambda > kappa
function result = loglikelihood_segment(data, lambda, kappa, beta)
	if (lambda <= kappa)
		result = -Inf;
		return;
	end
	X = data(1,:);
	Y = data(2,:);
	theta = -atan(beta);
	X2=X*cos(theta)-Y*sin(theta);
	minX=min(X2);
	maxX=max(X2);
	if (minX < kappa || maxX > lambda)
		result = -Inf;
		return;
	end
	result = -log(lambda-kappa);
end

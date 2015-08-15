function result = loglikelihood_line(data, alpha, beta, sigma)
	X = data(1,:);
	Y = data(2,:);
	Y_model = alpha + beta * X;
	N = size(X, 2);
	result = -N/2* (log(2*pi*sigma^2)) - sum((Y-Y_model).^2)/(2*sigma^2) ;
end


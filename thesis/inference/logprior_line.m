function result = logprior_line(alpha, beta, sigma)
	if (sigma < 0)
		result = -Inf;
		return;
	end
	result = -log(sigma) + log(1 + beta^2) * (-3/2);
end



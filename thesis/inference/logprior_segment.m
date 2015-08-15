% two-sided Pareto prior
% here "a" must be "outside of the prior" for two endpoints
% we use for the hyperparameters the greek symbols: xi, nu
% we assume xi > nu, lambda > kappa
function result = logprior_segment(lambda, kappa, xi, nu, K)
	if (xi <= nu)
		result = -Inf;
		return;
	end

	% shift to origin, now l and -l are endpoints, so only use +l
	p_shift = (xi + nu) / 2;
	xi_shift = xi - p_shift;
	lambda_shift = lambda - p_shift;
	kappa_shift = kappa - p_shift;

	if (lambda_shift >= xi_shift && kappa_shift <= -xi_shift)
		result = -(K+1) * (log(lambda_shift) + log(-kappa_shift)) + 2*K * log(xi_shift) + 2*log(K/2);
		%result = -(K+1) * log(lambda_shift) + K * log(xi_shift) + log(K/2) + -(K+1) * log(-kappa_shift) + K * log(xi_shift) + log(K/2);
	else
		result = -Inf;
	end
end

%logP2(11, -11, 10, -10, 3)



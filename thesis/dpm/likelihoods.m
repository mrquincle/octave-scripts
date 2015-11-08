function n = likelihoods(S, data, mu, Sigma)
	if (S.prior == 'NIW')
		n = exp(loggausspdf(data, mu, Sigma));
	elseif (S.prior == 'NIG')
	%	data
		y=data(end,:)';
		X=data(1:end-1,:);
		%X0=data(1:end-1,:);
		%X=[ones(1,size(X0,2)); X0];
	%	size(X)
	%	size(mu)
		mu2=sum(X.*mu',1)';
		n = exp(loggausspdf(y, mu2, Sigma));
	else
		error("Unknown type of prior");
	end
end

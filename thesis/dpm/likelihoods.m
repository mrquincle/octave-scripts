function n = likelihoods(S, data, mu, Sigma)
	switch(S.prior)
	case 'NIW'
		n = exp(loggausspdf(data, mu, Sigma));
	case 'NIG'
	%	data
		y=data(end,:)';
		X=data(1:end-1,:);
		%X0=data(1:end-1,:);
		%X=[ones(1,size(X0,2)); X0];
	%	size(X)
	%	size(mu)
		mu2=sum(X.*mu',1)';
		n = exp(loggausspdf(y, mu2, Sigma));
	otherwise
		error("Unknown type of prior");
	end
end

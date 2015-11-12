function R = sample_pdf(S)

	switch(S.prior)
	case 'NIW'
		[R.mu, R.Sigma] = normalinvwishrnd(S);
	case 'NIG'
		[R.mu, R.Sigma] = nigrnd(S);
	otherwise
		error("Unknown type of prior");
	end
end

% Sample probability density function
function R = sample_pdf(S)

	switch(S.prior)
	case 'NIW'
		[R.mu, R.Sigma] = normalinvwishrnd(S);
	case 'NIG'
		[R.mu, R.Sigma] = nigrnd(S);
	case 'DPM_Seg'
		[R.mu, R.Sigma] = nigrnd(S);
		[R.a, R.b] = pareto2rnd(S);
	otherwise
		error('Unknown type of prior');
	end
end

function R = sample_pdf(S)

	if (S.prior == 'NIW')
		[R.mu, R.Sigma] = normalinvwishrnd(S);
	elseif (S.prior == 'NIG')
		[R.mu, R.Sigma] = nigrnd(S);
	else
		error("Unknown type of prior");
	end
end

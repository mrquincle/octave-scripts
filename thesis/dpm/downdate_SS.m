function R = downdate_SS(z, S)

	if (S.prior == 'NIW')
		R = niwdowndate(z, S);
	elseif (S.prior == 'NIG')
		R = nigdowndate(z, S);
	else
		error("Unknown type of prior");
	end

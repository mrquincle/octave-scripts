function R = downdate_SS(z, S)
	switch (S.prior)
	case 'NIW'
		R = niwdowndate(z, S);
	case 'NIG'
		R = nigdowndate(z, S);
	otherwise
		error("Unknown type of prior");
	end
end

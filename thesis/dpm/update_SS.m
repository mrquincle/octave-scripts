function R = update_SS(z, S)
	% Update sufficient statistics with new data z
	% The data can be just observations itself or in the case of a regression problem it can be pairs of (X,y) with X the independent variables (can be a vector) and y the dependent variable

	if (S.prior == 'NIW')
		R = niwupdate(z, S);
	elseif (S.prior == 'NIG')
		R = nigupdate(z, S);
	else
		error("Unknown type of prior");
	end

end

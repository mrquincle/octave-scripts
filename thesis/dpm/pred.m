function out = pred(z, S)

% This is the posterior predictive

% it is NOT the marginal likelihood

	if (S.prior == 'NIW')

		% Get marginal likelihood for gaussian distribution with mean and
		% covariance matrix being distributed from a Normal inverse Wishart
		% distribution with parameters given by S

		mu0 = S.mu;
		kappa0 = S.kappa;
		nu0 = S.nu;
		lambda0 = S.lambda;

		p = length(mu0);

		% different S calculated here...
		%S = lambda0*(kappa0+1)/kappa0/(nu0-p+1);
		%nu = nu0-p+1;
		%out = (1+(z-mu0)'*S^-1*(z-mu0)/nu)^(-(nu+p)/2)*exp(gammaln((nu+p)/2))/...
		%    exp(gammaln((nu)/2))*(det(nu*p*S))^-.5;

		%n=p-1;
		%nu_n = nu0+n;
		%exp(gammaln(nu_n/2))/exp(gammaln(nu0/2))
		%(kappa_0/kappa_n)^(p/2);
		%lambda

		% looks like n=1

		S = lambda0*(kappa0+1)/kappa0/(nu0-p+1);
		nu = nu0-p+1;
		out = (1+(z-mu0)'*S^-1*(z-mu0)/nu)^(-(nu+p)/2)*...
			exp(gammaln((nu+p)/2))/exp(gammaln((nu)/2))*...
			(det(nu*p*S))^-.5;

		%snupi=det(nu*p*S)
		%c0=exp(gammaln((nu+p)/2))/exp(gammaln((nu)/2))
		%c=c0*(snupi)^-.5
		%diff=z-mu0
		%term=diff'*S^-1*diff
		%scatter=(1+term/nu)^(-(nu+p)/2)
		%result=scatter*c


	elseif (S.prior == 'NIG')
		% p(z*|Z) = Student_nu(X*mu, b/a(I + X Gamma^-1 X))
		% all z are summarized through sufficient statistics, so we need only
		% p(z*|ss) and here z* consists of y* and X*

		y=z(end,:);
		X=z(1:end-1,:);
		nu = 2*S.a;
		XGX=X' * S.Lambda * X;
		out = mvtpdf(y, X'*S.mu, S.b / S.a * (eye(size(XGX)) + XGX), nu);
	else
		error("Unknown type of prior");
	end

end

% Sample from a two-sided Pareto distribution
%
% -- Function: rnd_two_sided_pareto(x_l, x_m, alpha)
%     Return a sample from a two_sided Pareto distribution using inverse
%     transform sampling.
% -- Function: rnd_two_sided_pareto(x_l, x_m, alpha, R, C)
%     Return samples from a two-sided Pareto distribution in R rows and C
%     columns.
%
%     Scale parameters: x_l < 0, x_m > 0.
%     Shape parameter: alpha > 0.
%
%     The two-sided Pareto distribution can be used as a prior for a uniform
%     distribution over a segment.
%
%     The symmetric Pareto distribution can be sampled from by choosing
%     x_l = -x_m.

function sample = rnd_two_sided_pareto(x_l, x_m, alpha, R=1, C=1)
	posneg = binornd(1, 0.5, R, C);

	x_lm = posneg.*x_l + (1-posneg).*x_m;

	r = unifrnd(0, 1, R, C);
	sample = 1./(r.^(1/alpha)) .* x_lm;
end

% Sample from a Pareto distribution
%
% -- Function: rnd_pareto(x_m, alpha)
%     Return a sample from a Pareto distribution using inverse transform
%     sampling.
% -- Function: rnd_pareto(x_m, alpha, R, C)
%     Return samples from a Pareto distribution in R rows and C columns.
%
%     Scale parameter x_m is real, and larger than 0. For x < x_m, p(x) is
%     zero. Shape parameter alpha is real, and larger than 0. It defines the
%     probability p(x=x_m)=alpha.
%
%     Remark: unifrnd is used, according to help it says that it generates
%     samples from [A, B], but it seems to use rand() which draws from (A, B).
%     In our inverse transform sampling step we actually need to draw from
%     (A, B].
%
%     In inference for a uniform distribution with a Pareto prior, you will
%     have N, the number of samples, and x_m should not just be the extreme,
%     but max(x_m, prior).

function S = rnd_pareto(x_m, alpha, N=0, R=1, C=1)
	u = unifrnd(0, 1, R, C);
	S = u.^(-1/(alpha+N)) .* x_m;
end


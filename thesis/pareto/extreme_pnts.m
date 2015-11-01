% Get the most extreme pair of points from 1-dimensional or 2-dimensional
% points.
%
% -- Function: [A B] = extreme_points(Y)
%     Return a minimum point A and a maximum point B from the set of points Y.
%     The 'deterministic' method is used by default.
% -- Function: [A B] = extreme_points(Y, method)
%     Return a minimum point A and a maximum point B from the set of points Y.
%     The method can be 'deterministic' in which the first coordinate defines
%     what constitutes the minimum or maximum. Or 'coupled' in which any
%     coordinate can function as the one that defines the dimension over which
%     it is decided what constitutes the minimum or maximum. Or 'random' in
%     which this coordinate is randomly chosen per point.

function [A B] = extreme_pnts (Y, method='deterministic')
	dim = size(Y(1),2);

	if (dim == 1)
		A = min(Y);
		B = max(Y);
	else if (dim == 2)
		[w iw] = min(Y);
		e0 = Y(iw(1),:);
		e1 = Y(iw(2),:);
		[w iw] = max(Y);
		e2 = Y(iw(1),:);
		e3 = Y(iw(2),:);
		if (method == 'deterministic')
			A = e0;
			B = e2;
		elseif (method == 'coupled')
			warning("To implement. A and B need to be chosen at random");
			if (unidrnd(2)-1)
				A = e0;
				B = e2;
			else
				A = e1;
				B = e3;
			end
		elseif (method == 'random')
			if (unidrnd(2)-1)
				A = e0;
			else
				A = e1;
			end
			if (unidrnd(2)-1)
				B = e2;
			else
				B = e3;
			end
		end
	else
		warning("To implement. Get extreme points in high-dimensional spaces");
	end
end

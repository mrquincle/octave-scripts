% plots a line a+b*x from x0 to x1
function plot_line_segment(x0, x1, a, b)
	% angle, reverse
	theta = atan(b);
	%v = abs(a/b);
	if (b > 0)
		v = a/b;
	else
		v = -a/b;
	end
	p0 = [x0; v];
	p1 = [x1; v];
	R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
	D0 = R*p0;
	D1 = R*p1;
	x0=D0(1,:);
	y0=D0(2,:);
	x1=D1(1,:);
	y1=D1(2,:);
	plot([x0 x1], [y0 y1], 'k-');
	hold on
	plot(x0, y0, 'r*');
	plot(x1, y1, 'g*');
end

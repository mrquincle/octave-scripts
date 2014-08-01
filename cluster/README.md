
# Different scripts around the topic of clustering

## weights.m

Creating a series of points on a line and assign different weights to them. 

## lines.m

A more flexible octave file that allows you to specify a variable number of line coefficients. Subsequently points
according to these lines are drawn from the line equations itself, or with added white noise to both the x and the y 
coordinate depending on `config_use_noise`.

The script has several test functions:

    test_vertical_line: adjusts the second point such that is exactly above the first
    test_horizontal_line: adjusts the third point such that is exactly left from the first
    ...

To get from the points in `F` to all possible differences between the rows in `F` (with an `x` and `y` coordinate), the
function `Fdiff = bsxfun(@minus, F, permute(F, [3, 2, 1]))` is used. Nice to remember if you need it to operate on 
pairwise rows or columns some time. The angle is calculated as follows: `A = atan2(Fdiff(:,2,:), Fdiff(:,1,:))`. So,
this is the angle for the line between each pair of points. The angle between a point with itself is 0. The angle of a
point with a point vertical above it is `pi/2`. Due to the fact that we do not care about the "order" of the points, we
can map the angle of `pi*3/2` of a vertical line between a point and a point below it, to `pi/2` as well. Note, that 
with an angle we can present all possible lines with a value between `0` and `pi`, also vertical ones which have an
infinite slope.

Now, the distance to the origin has to be calculated for the line that is drawn to pairwise points. The distance of a
line with a point is always considered to be the distance of the closest point on a line with that point, in case you
are wondering.

If you check any source, e.g. [Wikipedia](https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line) you will see
that the equation for the distance to the origin is:

<!--
http://latex.codecogs.com/gif.latex?d=\frac{|-x_1y_2+x_2y_1|}{\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}}
-->

![equation](http://latex.codecogs.com/gif.latex?d%3D%5Cfrac%7B%7C-x_1y_2%2Bx_2y_1%7C%7D%7B%5Csqrt%7B%28x_2-x_1%29%5E2%2B%28y_2-y_1%29%5E2%7D%7D)

This is unfortunate because the denominator can become `0` if we happen to compare a point with itself! We do not 
care so much about the distance to the origin for these points, but they should blow up suddenly to infinity. It can
be solved easily by realizing that we can rotate the two points on our estimated line around the origin! By rotating
the points exactly by the estimated angle, we will have a horizontal line which will always pass through the `y`-axis.
Not only that, `|y|` is exactly the distance with the origin we are after. To see that you will have draw a few example
points.

This script does also some line estimation where you have to weigh the different points of each line yourself. It is 
illustrative to see how much influence a few pseudo-outliers (outliers in the form of points on another line) have on
the estimated line parameters. This line estimation is by the way orthogonal least squares, also known as Demming 
regression. And then in a weighted form. I can tell you beforehand, this sucks big time.

## deming.m

This is script is not by own. Check the implementation in `orthogonal.m` if you don't like his license.

## orthogonal.m

This script is a bit less advanced than the `lines.m` script. It only accepts parameters for one line. It shows four
figures of the successive steps of orthogonal least squares (or deming regression) when outliers are involved. In case
the outliers are assigned lower weight, the results are definitely reasonable. However, if there is only one outlier
with a weight as large as the normal points, the entire line estimation becomes very bad.




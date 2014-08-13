# Nonparametric Bayesian methods

To do inference in man-made environments which will be able to do inference in an environment such as a supermarket for our robots, we would like to be able to generate those structures first.

This is a fun part of generative Bayesian methods, which brings together nice mathematical models and intuition.

As an example of a model that generates lines according to a Hierarchical Dirichlet Process give and take a few tweaks you can check out [dpmobjectrnd.m](https://github.com/mrquincle/octave-scripts/blob/master/thesis/dpmobjectrnd.m). The lines are not all chosen randomly, but there is a Chinese Restaurant Process (CRP) which generates a finite number of orientations and lengths for each `line type`, hereby resembling something akin to a real life process where an object might generate parallel lines with a much greater chance than if there were no higher level objects present.

![manmade](pictures/manmade.png?raw=true "Man-made generation of lines")

The figure does of course not resemble something you recognize, but you see that this might be a nice starting point for inference if our inference method starts reasoning with this type structure in mind.

It is not so sophisticated to extends this to 2D parallelograms.

<p align="center">
<img src="pictures/patches.png?raw=true" alt="Man-made generation of patches" height="500px"/>
</p>

These patches do not resemble 3D wireframes, but you see that this is in the right direction.

Frames can be added by an additional random probability which picks on of the lines corresponding to the patch, this can be turned on by using `extremes_only` in the `dpmsquare1.config` file.

<p align="center">
<img src="pictures/frames.png?raw=true" alt="Man-made generation of frames" height="500px"/>
</p>

It is now clear as day that there are correlations between the lines, namely they are "glued together" in the form of squares. This is pretty hard to do if you do not start from the representation of a square.

<p align="center">
<img src="pictures/cubes.png?raw=true" alt="Man-made generation of cubes" height="500px"/>
</p>

Something resembling wireframes can be seen in this picture. There are no conditions that prevent cubes from intersection. This does not need to be an immediate problem, intersecting boxes might be seen as a composed object. However, in the end, we probably would like to be able to formulate conditions such as boxes are likely to be stacked on top of each other, and other macro-object constraints that have to do with gravity and volumetric exclusion.

<p align="center">
<img src="pictures/regular4.png?raw=true" alt="Show regular patterns" height="500px"/>
</p>

The incorporation of regularity is quite hard. The figure above is one of the first attempts, but this needs much more attention. There are intricate things going on between the size of the objects, the spacing of a grid, and the number of objects in a grid. Generating for example the spacing as well as the number of items along one of the dimensions of grid from one multinomial distribution leads to on average grid sizes that have too many items, say 40x100, and sizes of 30x40 "pixels" (so, also a lot of overlap).

## Inference

To perform inference several methods can be used. One that fits especially discrete distribution is slice sampling. The following, however, is obtained through a collapsed CRP. First we show the lines we have randomly generated:

<p align="center">
<img src="pictures/lines_rnd.png?raw=true" alt="Generated lines" height="500px"/>
</p>

And here is the result from the inference process:

<p align="center">
<img src="pictures/lines_inf0.png?raw=true" alt="Inference of lines" height="500px"/>
</p>

The graphics are different because matlab is used now instead of octave. The algorithms from Caron are superslow in octave. Later on I will implement a fast inference method in C++, so this will not be the case in the near future. 

The inference procedure uses the same configuration file as the generation procedure for now. So, we do not search for the hyperparameter alpha for example. This would be normally part of the inference process. Note that a generation process with a high alpha will generate many short lines as well. Therefore, for the inference process to only capture the most "salient lines" it is wise to set alpha much lower than the "true" value. This is an interesting observation in my opinion. 

## Copyrights

* Author: Anne van Rossum
* License: LGPLv3
* Copyright: Distributed Organisms B.V. ([DoBots](http://dobots.nl)) and [Almende B.V.](http://almende.com)

Please, if you want to work on this problem, contact me. I would be sad if you just use it and I don't get published because someone else used my code and published just before me... So, I would really like to collaborate in that case. 



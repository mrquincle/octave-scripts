# Nonparametric Bayesian methods

To do inference in man-made environments which will be able to do inference in an environment such as a supermarket for our robots, we would like to be able to generate those structures first.

This is a fun part of generative Bayesian methods, which brings together nice mathematical models and intuition.

As an example of a model that generates lines according to a Hierarchical Dirichlet Process give and take a few tweaks you can check out [dpmobjectrnd.m](https://github.com/mrquincle/octave-scripts/blob/master/thesis/dpmobjectrnd.m). The lines are not all chosen randomly, but there is a Chinese Restaurant Process which generates a finite number of orientations and lengths for each `line type`, hereby resembling something akin to a real life process where an object might generate parallel lines with a much greater chance than if there were no higher level objects present.

![manmade](pictures/manmade.png?raw=true "Man-made generation of lines")

The figure does of course not resemble something you recognize, but you see that this might be a nice starting point for inference if our inference method starts reasoning with this type structure in mind.

It is not so sophisticated to extends this to 2D parallelograms.

![patches](pictures/patches.png?raw=true "Man-made generation of patches")

These patches do not resemble 3D wireframes, but you see that this is in the right direction.

## Inference

Currently, no inference is done. However, what I will do first is simple Gibbs sampling. After that I will develop an inference method that makes use of the hierarchical structure in the model.

## Copyrights

* Author: Anne van Rossum
* License: LGPLv3
* Copyright: Distributed Organisms B.V. ([DoBots](http://dobots.nl)) and [Almende B.V.](http://almende.com)

Please, if you want to work on this problem, contact me. I would be sad if you just use it and I don't get published because someone else used my code and published just before me... So, I would really like to collaborate in that case. 



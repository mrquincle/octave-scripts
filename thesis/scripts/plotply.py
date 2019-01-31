#!/usr/bin/python3

import numpy as np
import matplotlib.pylab as plt

from argparse import ArgumentParser

from mpl_toolkits.mplot3d import Axes3D
from third_party.python_plyfile.plyfile import PlyElement, PlyData

def load_ply(file_name):
    ply_data = PlyData.read(file_name)
    points = ply_data['vertex']
    print("File contains # points:", points.count)
    #points = np.vstack([points['x'], points['y'], points['z']]).T
    ret_val = [points]

    if len(ret_val) == 1:  # Unwrap the list
        ret_val = ret_val[0]

    return ret_val

def plot_3d_point_cloud(x, y, z, show=True, show_axis=True, in_u_sphere=False, marker='.', s=8, alpha=.8, figsize=(5, 5), elev=10, azim=240, axis=None, title=None, *args, **kwargs):

    if axis is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')        
    else:
        ax = axis
        fig = axis

    if title is not None:
        plt.title(title)

    sc = ax.scatter(x, y, z, marker=marker, s=s, alpha=alpha, *args, **kwargs)
    ax.view_init(elev=elev, azim=azim)

    if in_u_sphere:
        ax.set_xlim3d(-0.5, 0.5)
        ax.set_ylim3d(-0.5, 0.5)
        ax.set_zlim3d(-0.5, 0.5)
    else:
        miv = 0.7 * np.min([np.min(x), np.min(y), np.min(z)])  # Multiply with 0.7 to squeeze free-space.
        mav = 0.7 * np.max([np.max(x), np.max(y), np.max(z)])
        ax.set_xlim(miv, mav)
        ax.set_ylim(miv, mav)
        ax.set_zlim(miv, mav)
        plt.tight_layout()

    if not show_axis:
        plt.axis('off')

    if 'c' in kwargs:
        plt.colorbar(sc)

    if show:
        plt.show()

    return fig

def main():
    print("Plot .ply file")
    parser = ArgumentParser()
    parser.add_argument('ply_filename')

    args = parser.parse_args()
    print("Load file: ", args.ply_filename)
    points = load_ply(args.ply_filename)

    plot_3d_point_cloud(points['x'], points['y'], points['z'])

main()

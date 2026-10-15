#!/usr/bin/env python3
"""
check_stl_overlap.py

Check whether two STL meshes overlap and report:

    - overlap volume
    - maximum penetration depth of mesh A into mesh B

Usage:
    python check_stl_overlap.py mesh1.stl mesh2.stl

Optional:
    --samples N
        Number of surface samples used for penetration-depth calculation.
        Default: 100000

    --engine ENGINE
        Boolean engine for overlap-volume calculation.
        Examples: manifold, blender, scad

Notes:
    - Penetration depth is measured from the surface of mesh A into mesh B.
    - A positive penetration depth means that part of A's surface lies inside B.
    - The reported depth is approximate because the A surface is sampled.
    - Mesh B should be watertight for signed-distance calculations.
"""

import argparse
import sys
import os
import numpy as np
import trimesh


def load_mesh(path):
    mesh = trimesh.load_mesh(
        path,
        force="mesh",
        process=True
    )

    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError(
            "{} did not load as a single mesh object".format(path)
        )

    return mesh


def aabb_overlap(mesh_a, mesh_b):
    """Fast bounding-box rejection test."""
    min_a, max_a = mesh_a.bounds
    min_b, max_b = mesh_b.bounds

    return bool(
        np.all(max_a >= min_b) and
        np.all(max_b >= min_a)
    )


def combine_result(result):
    """Convert a boolean result into a single Trimesh."""
    if result is None:
        return None

    if isinstance(result, trimesh.Trimesh):
        return result

    if isinstance(result, trimesh.Scene):
        meshes = [
            geom
            for geom in result.geometry.values()
            if isinstance(geom, trimesh.Trimesh)
        ]

        if not meshes:
            return None

        return trimesh.util.concatenate(meshes)

    if isinstance(result, (list, tuple)):
        meshes = [
            geom
            for geom in result
            if isinstance(geom, trimesh.Trimesh)
        ]

        if not meshes:
            return None

        return trimesh.util.concatenate(meshes)

    raise TypeError(
        "Unsupported boolean result type: {}".format(type(result))
    )


def intersection_volume(mesh_a, mesh_b, engine=None):
    """
    Calculate the intersection volume.
    """
    result = trimesh.boolean.intersection(
        [mesh_a, mesh_b],
        engine=engine
    )

    result = combine_result(result)

    if result is None or len(result.faces) == 0:
        return 0.0

    volume = abs(float(result.volume))

    return volume


def maximum_penetration_depth(mesh_a, mesh_b, num_samples=100000):
    """
    Estimate the maximum depth that the surface of mesh A penetrates
    into mesh B.

    Positive depth means a sampled point on A is inside B.

    The calculation uses:
        signed_distance(B, points)

    where positive values are inside B and negative values are outside
    according to trimesh's signed-distance convention. To avoid depending
    on sign conventions, we explicitly use B.contains() to identify
    interior points and then measure their distance to B's surface.
    """

    if not mesh_b.is_watertight:
        raise ValueError(
            "Mesh B must be watertight to calculate penetration depth."
        )

    # Start with all vertices of A.
    points = mesh_a.vertices

    # Add surface samples if requested.
    if num_samples > len(points):
        sample_count = num_samples - len(points)

        sampled_points, _ = trimesh.sample.sample_surface(
            mesh_a,
            sample_count
        )

        points = np.vstack([
            points,
            sampled_points
        ])

    # Find which sampled points are actually inside B.
    inside = mesh_b.contains(points)

    if not np.any(inside):
        return 0.0

    inside_points = points[inside]

    # Distance from points inside B to the surface of B.
    closest_points, distances, _ = trimesh.proximity.closest_point(
        mesh_b,
        inside_points
    )

    if len(distances) == 0:
        return 0.0

    return float(np.max(distances))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Check STL mesh overlap, calculate overlap volume, "
            "and estimate maximum penetration depth for all .stl files in folder"
        )
    )

    parser.add_argument(
        "folder",
        help="Folder of .stl files"
    )

    parser.add_argument(
        "--samples",
        type=int,
        default=100000,
        help=(
            "Number of surface samples used for penetration "
            "calculation. Default: 100000"
        )
    )

    parser.add_argument(
        "--engine",
        default=None,
        help=(
            "Boolean engine for volume calculation: "
            "manifold, blender, or scad. Default: auto."
        )
    )

    args = parser.parse_args()

    if args.samples < 0:
        parser.error("--samples must be >= 0")

    for file1 in os.listdir(args.folder):
        if not file1.endswith('.stl'):
            continue    
        mesh_a = load_mesh(args.folder+"/"+file1)
        if(mesh_a.is_watertight):
            print(file1)  
        else:
            print(file1+" NOT WATERTIGHT")  
        for file2 in os.listdir(args.folder):
            if (not file1.endswith('.stl')) or file1==file2:
                continue           
            mesh_b = load_mesh(args.folder+"/"+file2)
            if not aabb_overlap(mesh_a, mesh_b):
                d=0         
#                print("    "+file2+"    NO")
            else:
                print("    OVERLAPS "+file2+" by volume: "+str(intersection_volume(mesh_a,mesh_b,args.engine))+
                                                                                                  ", max depth: "+str(maximum_penetration_depth(mesh_a,mesh_b)))

if __name__ == "__main__":
    main()

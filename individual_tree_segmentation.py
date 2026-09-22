# ============================================================
# Li2012 Individual Tree Segmentation with Plot Auto-Selection
# ============================================================

import os
import laspy
import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.geometry import MultiPoint
from shapely.ops import unary_union

from shapely.geometry import Point
from shapely import contains_xy

from scipy.spatial import cKDTree
from tqdm import tqdm



# ============================================================
# PARAMETERS
# ============================================================

INPUT_LAS = r"01_data/site 3.laz"

PLOTS_FILE = r"01_data/Baringo_Sites.geojson"

CRS_EPSG = 21097


OUTPUT_DIR = "tree_segmentation_output/site3"


# Li et al. 2012 parameters

DT1 = 2.0
DT2 = 1.5
ZU = 15


MIN_TREE_HEIGHT = 2

MIN_POINTS_TREE = 20



# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


TREE_DIR = os.path.join(
    OUTPUT_DIR,
    "individual_trees"
)


os.makedirs(
    TREE_DIR,
    exist_ok=True
)



# ============================================================
# LOAD PLOTS
# ============================================================

print("Loading plot polygons...")


plots = gpd.read_file(
    PLOTS_FILE
)



# create plot id

if "plot_id" not in plots.columns:

    plots["plot_id"] = np.arange(
        1,
        len(plots)+1
    )



# CRS correction

if plots.crs is None:

    plots = plots.set_crs(
        epsg=CRS_EPSG
    )


elif plots.crs.to_epsg() != CRS_EPSG:

    plots = plots.to_crs(
        epsg=CRS_EPSG
    )



print(
    "Number of plots:",
    len(plots)
)



# ============================================================
# READ LAS
# ============================================================

print("Reading LAS...")


las = laspy.read(
    INPUT_LAS
)



xyz = np.vstack(
    (
        las.x,
        las.y,
        las.z
    )
).T



print(
    "Original points:",
    len(xyz)
)



# ============================================================
# SELECT CORRESPONDING PLOT
# (Equivalent to choose_plot_by_centroid() in R)
# ============================================================


def choose_plot_by_centroid(
        xyz,
        plots
):


    xmin = xyz[:,0].min()
    xmax = xyz[:,0].max()

    ymin = xyz[:,1].min()
    ymax = xyz[:,1].max()



    centroid = Point(
        (xmin+xmax)/2,
        (ymin+ymax)/2
    )



    inside = plots.contains(
        centroid
    )


    hits = plots[
        inside
    ]



    if len(hits) > 0:

        print(
            "Centroid falls inside plot:",
            hits.iloc[0].plot_id
        )


        return hits.iloc[0].geometry



    # otherwise nearest plot

    distances = plots.geometry.distance(
        centroid
    )


    idx = distances.idxmin()


    print(
        "Nearest plot selected:",
        plots.loc[idx,"plot_id"]
    )


    return plots.loc[
        idx,
        "geometry"
    ]



polygon = choose_plot_by_centroid(
    xyz,
    plots
)



# ============================================================
# CLIP POINT CLOUD TO PLOT
# ============================================================


print("Cropping point cloud...")


xmin,ymin,xmax,ymax = polygon.bounds



# First bounding box crop

bbox_mask = (

    (xyz[:,0] >= xmin) &
    (xyz[:,0] <= xmax) &
    (xyz[:,1] >= ymin) &
    (xyz[:,1] <= ymax)

)



xyz_crop = xyz[
    bbox_mask
]



print(
    "After bbox crop:",
    len(xyz_crop)
)



# Exact polygon crop

inside = contains_xy(
    polygon,
    xyz_crop[:,0],
    xyz_crop[:,1]
)



xyz_crop = xyz_crop[
    inside
]



print(
    "After polygon crop:",
    len(xyz_crop)
)



if len(xyz_crop)==0:

    raise ValueError(
        "No points inside selected plot. Check CRS."
    )



# ============================================================
# HEIGHT NORMALIZATION
# ============================================================


print("Normalizing heights...")


if hasattr(
    las,
    "classification"
):


    classifications = np.array(
        las.classification
    )


    # apply same cropping mask

    class_crop = classifications[
        bbox_mask
    ][inside]



    ground = xyz_crop[
        class_crop == 2
    ]



    if len(ground) > 0:


        ground_tree = cKDTree(
            ground[:,:2]
        )


        _,idx = ground_tree.query(
            xyz_crop[:,:2],
            k=1
        )


        ground_z = ground[
            idx,
            2
        ]


        height = (
            xyz_crop[:,2]
            -
            ground_z
        )



    else:

        print(
            "No ground points found, using absolute height"
        )

        height = xyz_crop[:,2]



else:

    print(
        "No classification found"
    )

    height = xyz_crop[:,2]




# Remove low vegetation

mask = (
    height >= MIN_TREE_HEIGHT
)


xyz_crop = xyz_crop[
    mask
]


height = height[
    mask
]



print(
    "Vegetation points:",
    len(xyz_crop)
)



# ============================================================
# BUILD KD TREE
# ============================================================


kdtree = cKDTree(
    xyz_crop[:,:2]
)



# ============================================================
# TREE TOP DETECTION
# ============================================================


def detect_tree_tops(
        xyz,
        height,
        radius
):

    tops=[]


    for i,p in enumerate(
        tqdm(
            xyz,
            desc="Detecting tree tops"
        )
    ):


        neighbours = kdtree.query_ball_point(
            p[:2],
            radius
        )


        if height[i] >= np.max(
            height[neighbours]
        ):

            tops.append(i)



    return np.array(tops)



tree_tops = detect_tree_tops(
    xyz_crop,
    height,
    DT1
)



print(
    "Tree tops detected:",
    len(tree_tops)
)



tree_tops = tree_tops[
    np.argsort(
        height[tree_tops]
    )[::-1]
]



# ============================================================
# LI2012 REGION GROWING
# ============================================================


print("Running Li2012 segmentation...")


tree_ID = np.zeros(
    len(xyz_crop),
    dtype=np.uint32
)


assigned=set()

tree_counter=1



for seed in tqdm(
    tree_tops,
    desc="Growing crowns"
):


    if seed in assigned:
        continue


    queue=[seed]

    crown=[]



    while queue:


        current = queue.pop()



        if current in assigned:
            continue



        assigned.add(
            current
        )


        crown.append(
            current
        )



        radius = DT2


        if height[current] > ZU:

            radius = DT1



        neighbours = kdtree.query_ball_point(
            xyz_crop[current,:2],
            radius
        )



        for n in neighbours:


            if n in assigned:
                continue



            dz = abs(
                height[current]
                -
                height[n]
            )


            if dz < 3:

                queue.append(n)



    if len(crown) >= MIN_POINTS_TREE:


        tree_ID[crown] = tree_counter

        tree_counter += 1



print(
    "Trees segmented:",
    tree_counter-1
)



# ============================================================
# SAVE SEGMENTED CLOUD
# ============================================================


print("Saving segmented cloud...")


out_las = laspy.create(
    point_format=las.header.point_format,
    file_version=las.header.version
)


out_las.x = xyz_crop[:,0]
out_las.y = xyz_crop[:,1]
out_las.z = xyz_crop[:,2]


out_las.add_extra_dim(
    laspy.ExtraBytesParams(
        name="TreeID",
        type=np.uint32
    )
)


out_las.TreeID = tree_ID



out_las.write(
    os.path.join(
        OUTPUT_DIR,
        "segmented_plot.laz"
    )
)
# ============================================================
# CREATE TREE CROWN POLYGONS
# ============================================================


def create_tree_crowns(
        xyz,
        tree_ID,
        height
):

    polygons = []
    attributes = []


    unique_trees = np.unique(
        tree_ID[tree_ID > 0]
    )


    print(
        "Creating crown polygons..."
    )


    for tid in tqdm(
        unique_trees,
        desc="Crown polygons"
    ):


        idx = np.where(
            tree_ID == tid
        )[0]


        points = xyz[idx,:2]


        if len(points) < 3:
            continue



        multipoint = MultiPoint(
            points
        )


        # convex hull (robust)
        crown = multipoint.convex_hull



        if crown.is_empty:
            continue



        attributes.append(

            {
                "TreeID": int(tid),

                "Height":
                float(
                    np.max(height[idx])
                ),

                "Points":
                int(
                    len(idx)
                ),

                "CrownArea":
                float(
                    crown.area
                ),

                "CrownWidth":
                float(
                    np.sqrt(
                        crown.area
                    )
                ),

                "X":
                float(
                    np.mean(points[:,0])
                ),

                "Y":
                float(
                    np.mean(points[:,1])
                )
            }

        )


        polygons.append(
            crown
        )



    return polygons, attributes


# ============================================================
# EXPORT INDIVIDUAL TREES
# ============================================================


metrics=[]


print("Exporting trees...")


for tid in tqdm(
    np.unique(tree_ID[tree_ID>0])
):


    idx = np.where(
        tree_ID==tid
    )[0]


    pts = xyz_crop[idx]



    tree_las = laspy.create(
        point_format=las.header.point_format,
        file_version=las.header.version
    )


    tree_las.x = pts[:,0]
    tree_las.y = pts[:,1]
    tree_las.z = pts[:,2]



    tree_las.write(
        os.path.join(
            TREE_DIR,
            f"tree_{tid:04d}.laz"
        )
    )



    metrics.append(

        {
            "TreeID":tid,
            "Height":float(
                np.max(height[idx])
            ),
            "Points":len(idx),
            "X":float(
                np.mean(pts[:,0])
            ),
            "Y":float(
                np.mean(pts[:,1])
            )
        }

    )



# ============================================================
# SAVE METRICS
# ============================================================


pd.DataFrame(
    metrics
).to_csv(

    os.path.join(
        OUTPUT_DIR,
        "tree_metrics.csv"
    ),

    index=False

)
# ============================================================
# EXPORT TREE CROWN GEOPACKAGE
# ============================================================


crowns, crown_attributes = create_tree_crowns(
    xyz_crop,
    tree_ID,
    height
)



tree_gdf = gpd.GeoDataFrame(
    crown_attributes,
    geometry=crowns,
    crs=f"EPSG:{CRS_EPSG}"
)



gpkg_file = os.path.join(
    OUTPUT_DIR,
    "tree_segmentation.gpkg"
)



tree_gdf.to_file(
    gpkg_file,
    layer="tree_crowns",
    driver="GPKG"
)



print(
    "GeoPackage written:",
    gpkg_file
)


print("===================================")
print("TREE SEGMENTATION COMPLETE")
print(
    "Trees:",
    tree_counter-1
)
print("===================================")
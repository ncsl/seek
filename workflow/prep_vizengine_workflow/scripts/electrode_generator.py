import os
import os.path
import random
import sys

import bpy
import numpy as np
import pandas as pd


def text_num_split(item):
    """Split string into alphabet and numeric parts.

    Requires that string starts with alphabet and ends
    with numeric char.
    """
    for index, letter in enumerate(item, 0):
        if letter.isdigit():
            return [item[:index], item[index:]]


def main(elec_fpath):
    fs_subjects_dir = os.environ.get("SUBJECTS_DIR")
    subject = os.environ.get("SUBJECT")
    subjects_dir = f"{fs_subjects_dir}/{subject}"

    scn = bpy.context.scene
    if not scn.render.engine == "CYCLES":
        scn.render.engine = "CYCLES"

    # create blender object
    bpy.ops.object.empty_add()
    bpy.context.active_object.name = "Electrodes"

    elec_df = pd.read_csv(elec_fpath, sep="\t")

    # keep track of which electrode group
    oldElectrodeGroup = ""

    # loop through all electrodes stored in electrodes.tsv file
    for index, elecName in enumerate(elec_df["name"]):
        # if left hemisphere, then there will be a "'" character
        # else, right hemisphere
        print(f'Creating blender object for {elecName}...')
        if "'" in elecName:
            electrodeGroup = elecName.split("'")[0]
        else:
            electrodeGroup = text_num_split(elecName)[0]

        # create a new electrode group
        if oldElectrodeGroup != electrodeGroup:
            bpy.ops.object.empty_add()
            bpy.context.active_object.name = electrodeGroup
            bpy.context.active_object.parent = bpy.data.objects["Electrodes"]
            oldElectrodeGroup = electrodeGroup
            mat = bpy.data.materials.new("electrodeMaterial")
            mat.diffuse_color = (
                random.random(),
                random.random(),
                random.random(),
                1,
            )
            mat.use_nodes = True
        # electrodeName = "{group}{name}".format(
        #     group=electrodeGroup, name=elecName.split("'")[1]
        # )

        # check of nans / 'n/a' in elec_df
        if any(np.isnan(elec_df[col][index]) for col in ["x", "y", "z"]):
            continue

        electrodeX = float(elec_df["x"][index])
        electrodeY = float(elec_df["y"][index])
        electrodeZ = float(elec_df["z"][index])

        # create a blender-mesh sphere at the electrode location
        bpy.ops.mesh.primitive_ico_sphere_add(
            location=(electrodeX, electrodeY, electrodeZ)
        )
        bpy.context.active_object.name = elecName
        bpy.context.active_object.active_material = mat
        bpy.context.active_object.parent = bpy.data.objects[electrodeGroup]

    bpy.ops.export_scene.gltf(
        export_format="GLB",
        filepath=f"{subjects_dir}/blender_objects/electrodes",
        export_texcoords=False,
        export_normals=False,
        export_cameras=False,
        export_yup=False,
    )


if __name__ == "__main__":
    elec_fpath = sys.argv[6]

    main(elec_fpath)

import sys
import json
import os

import bpy


def main(mat_color_json_fpath):
    patientDir = os.environ.get("SUBJECTS_DIR")
    patientID = os.environ.get("SUBJECT")

    print(f'Loading materialColors.json file from {mat_color_json_fpath}')
    with open(mat_color_json_fpath, 'r', encoding='utf-8') as json_file:
        data = json.load(json_file)

    scn = bpy.context.scene
    if not scn.render.engine == "CYCLES":
        scn.render.engine = "CYCLES"

    subjDir = "{patientDir}/{patient}".format(patientDir=patientDir, patient=patientID)

    bpy.ops.object.empty_add()
    bpy.context.active_object.name = "Brain"
    bpy.ops.object.empty_add()
    bpy.context.active_object.name = "Gyri"
    bpy.context.active_object.parent = bpy.data.objects["Brain"]
    bpy.ops.object.empty_add()
    bpy.context.active_object.name = "WhiteMatter"
    bpy.context.active_object.parent = bpy.data.objects["Brain"]
    bpy.ops.object.empty_add()
    bpy.context.active_object.name = "SubcorticalStructs"
    bpy.context.active_object.parent = bpy.data.objects["Brain"]

    for file in os.listdir("{dir}/obj".format(dir=subjDir)):
        name = file.split(sep=".obj")[0]

        r = float(data[os.path.splitext(file)[0]][0] / 255)
        g = float(data[os.path.splitext(file)[0]][1] / 255)
        b = float(data[os.path.splitext(file)[0]][2] / 255)
        if name.endswith("-Cerebral-Cortex") == False:
            bpy.ops.import_scene.obj(
                filepath="{dir}/obj/".format(dir=subjDir) + file, axis_forward="Y"
            )
            mat = bpy.data.materials.new("brainMaterial")
            mat.diffuse_color = (r, g, b, 1)
            bpy.data.objects[os.path.splitext(file)[0]].active_material = mat
            mat.use_nodes = True

            if name.endswith("-Cerebral-Cortex"):
                pass
            elif file[2] == ".":
                bpy.data.objects[os.path.splitext(file)[0]].parent = bpy.data.objects[
                    "Gyri"
                ]
            elif (
                file == "Left-Cerebral-White-Matter.obj"
                or file == "Right-Cerebral-White-Matter.obj"
            ):
                bpy.data.objects[os.path.splitext(file)[0]].parent = bpy.data.objects[
                    "WhiteMatter"
                ]
            else:
                bpy.data.objects[os.path.splitext(file)[0]].parent = bpy.data.objects[
                    "SubcorticalStructs"
                ]

    bpy.ops.export_scene.gltf(
        export_format="GLB",
        filepath="{dir}/{patient}".format(dir=subjDir, patient="brain"),
        export_texcoords=False,
        export_normals=False,
        export_cameras=False,
        export_yup=False,
    )


if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument(
    #     "--materialColorsPath",
    #     # dest="describe",
    #     # action="store_true",
    #     default="/home/scripts/materialColors.json",
    #     help=("If set, is the path to the materialColors JSON file."),
    # )
    #
    # args = parser.parse_args()

    # Extract arguments from parser
    mat_colors_json_fpath = sys.argv[5]
    if mat_colors_json_fpath == '':
        mat_colors_json_fpath = "/home/scripts/materialColors.json"

    main(mat_color_json_fpath=mat_colors_json_fpath)

import argparse
import json
import math
import os
import sys
import bpy


def main(
    fsdir,
    patient,
    electrode_file,
    fbx_output_fpath,
    glb_output_fpath,
    materialcolors_fpath,
    electrodeExport=False,
    justCortex=False,
):
    # print(fsdir, patient, electrode_file, fbx_output_fpath, glb_output_fpath, materialcolors_fpath)
    # print("trying to read: ", materialcolors_fpath)
    with open(materialcolors_fpath) as json_file:
        data = json.load(json_file)

    scn = bpy.context.scene
    if not scn.render.engine == "CYCLES":
        scn.render.engine = "CYCLES"

    subjDir = os.path.join(fsdir, "{patient}".format(patient=patient))
    if electrodeExport:
        bpy.ops.object.empty_add(type="CUBE")
        bpy.context.active_object.name = "Electrodes"
        bpy.context.active_object.rotation_euler = (-math.pi / 2, 0, math.pi)
        bpy.context.active_object.location = (128, 128, 128)
        electrodes = open(electrode_file)
        elecs = electrodes.readlines()
        for elec in elecs:
            electrodeGroup = elec.split("\t")[0]
            electrodeName = elec.split("\t")[2]
            # electrodeName = electrodeGroup + '_'+elec.split('\t')[2]
            electrodeX = float(elec.split("\t")[3])
            electrodeY = float(elec.split("\t")[4])
            electrodeZ = float(elec.split("\t")[5])
            bpy.ops.mesh.primitive_ico_sphere_add(
                location=(electrodeX, electrodeY, electrodeZ)
            )
            bpy.context.active_object.name = electrodeName
            bpy.context.active_object.parent = bpy.data.objects["Electrodes"]
    else:
        pass

    if justCortex != True:
        bpy.ops.object.empty_add(type="CUBE")
        bpy.context.active_object.name = "Brain"
        bpy.ops.object.empty_add(type="CUBE")
        bpy.context.active_object.name = "Gyri"
        bpy.ops.object.empty_add(type="CUBE")
        bpy.context.active_object.name = "WhiteMatter"
        bpy.context.active_object.rotation_euler = (0, 0, 0)

        for file in os.listdir("{dir}/obj".format(dir=subjDir)):
            if (
                file == "Right-Cerebral-Cortex.obj"
                or file == "Left-Cerebral-Cortex.obj"
            ):
                pass
            else:
                r = float(data[os.path.splitext(file)[0]][0] / 255)
                g = float(data[os.path.splitext(file)[0]][1] / 255)
                b = float(data[os.path.splitext(file)[0]][2] / 255)
                bpy.ops.import_scene.obj(
                    filepath="{dir}/obj/".format(dir=subjDir) + file
                )
                mat = bpy.data.materials.new("brainMaterial")
                mat.diffuse_color = (r, g, b, 1)
                bpy.data.objects[os.path.splitext(file)[0]].active_material = mat
                mat.use_nodes = True
                if file[2] == ".":
                    bpy.data.objects[
                        os.path.splitext(file)[0]
                    ].parent = bpy.data.objects["Gyri"]
                elif (
                    file == "Left-Cerebral-White-Matter.obj"
                    or file == "Right-Cerebral-White-Matter.obj"
                ):
                    bpy.data.objects[
                        os.path.splitext(file)[0]
                    ].parent = bpy.data.objects["WhiteMatter"]
                else:
                    bpy.data.objects[
                        os.path.splitext(file)[0]
                    ].parent = bpy.data.objects["Brain"]
                # bpy.data.objects['WhiteMatter'].parent = bpy.data.objects['Brain']
                # bpy.data.objects['Gyri'].parent = bpy.data.objects['Brain']
    else:
        bpy.ops.import_scene.obj(
            filepath="{dir}/obj/".format(dir=subjDir) + "Left-Cerebral-Cortex.obj"
        )
        mat = bpy.data.materials.new("brainMaterial")
        mat.diffuse_color = (float(245 / 255), float(245 / 255), float(245 / 255), 1)
        o = bpy.context.selected_objects[0]
        o.active_material = mat

        bpy.ops.import_scene.obj(
            filepath="{dir}/obj/".format(dir=subjDir) + "Right-Cerebral-Cortex.obj"
        )
        o = bpy.context.selected_objects[0]
        o.active_material = mat

    bpy.ops.export_scene.fbx(filepath=fbx_output_fpath)
    bpy.ops.export_scene.gltf(
        filepath=glb_output_fpath, export_texcoords=False, export_normals=False
    )


if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("fsdir", help="The freesurfer output diretroy.")
    # parser.add_argument("patid")
    # parser.add_argument(
    #     "electrode_initialization_file",
    #     help="The electrode file with contacts localized.",
    # )
    # parser.add_argument(
    #     "initial_run",
    #     help="Is this the first time running for patient?",
    # )
    # parser.add_argument(
    #     "just_cortex",
    #     help="Just do the cortex, or also include subcortex?",
    # )
    # parser.add_argument(
    #     "fbx_output_fpath",
    # )
    # parser.add_argument(
    #     "glb_output_fpath",
    # )
    # args = parser.parse_args()
    #
    # # extract arguments from parser
    # fsdir = args.fsdir
    # patid = args.patid
    # electrode_file = args.electrode_initialization_file
    # initial_run = args.initial_run
    # just_cortex = args.just_cortex
    # fbx_output_fpath = args.fbx_output_fpath
    # glb_output_fpath = args.glb_output_fpath

    fsdir = sys.argv[7]
    patid = sys.argv[8]
    electrode_file = sys.argv[9]
    initial_run = sys.argv[10]
    just_cortex = sys.argv[11]
    fbx_output_fpath = sys.argv[12]
    glb_output_fpath = sys.argv[13]
    materialcolors_fpath = sys.argv[14]

    print(sys.argv[:])
    main(
        fsdir,
        patid,
        electrode_file,
        fbx_output_fpath,
        glb_output_fpath,
        materialcolors_fpath,
        initial_run,
        just_cortex,
    )

import bpy
import os
import sys
import json
import math


def main(patient='', electrodeExport=False, justCortex=False):

    with open('scripts/materialColors.json') as json_file:
        data = json.load(json_file)

    scn = bpy.context.scene
    if not scn.render.engine == 'CYCLES':
        scn.render.engine = 'CYCLES'

    subjDir = "{dir}{patient}".format(
        dir=os.environ['SUBJECTS_DIR'], patient=patient)

    if electrodeExport:
        bpy.ops.object.empty_add(type='CUBE')
        bpy.context.active_object.name = 'Electrodes'
        bpy.context.active_object.rotation_euler = (-math.pi/2, 0, math.pi)
        bpy.context.active_object.location = (128, 128, 128)
        electrodes = open(
            "{dir}/electrodes/electrodes.txt".format(dir=subjDir))
        elecs = electrodes.readlines()
        for elec in elecs:
            electrodeGroup = elec.split('\t')[0]
            electrodeName = elec.split('\t')[2]
            # electrodeName = electrodeGroup + '_'+elec.split('\t')[2]
            electrodeX = float(elec.split('\t')[3])
            electrodeY = float(elec.split('\t')[4])
            electrodeZ = float(elec.split('\t')[5])
            bpy.ops.mesh.primitive_ico_sphere_add(
                location=(electrodeX, electrodeY, electrodeZ))
            bpy.context.active_object.name = electrodeName
            bpy.context.active_object.parent = bpy.data.objects['Electrodes']
    else:
        pass

    if justCortex != True:
        bpy.ops.object.empty_add(type='CUBE')
        bpy.context.active_object.name = 'Brain'
        bpy.ops.object.empty_add(type='CUBE')
        bpy.context.active_object.name = 'Gyri'
        bpy.ops.object.empty_add(type='CUBE')
        bpy.context.active_object.name = 'WhiteMatter'
        bpy.context.active_object.rotation_euler = (0, 0, 0)

        for file in os.listdir("{dir}/obj".format(dir=subjDir)):
            if file == 'Right-Cerebral-Cortex.obj' or file == 'Left-Cerebral-Cortex.obj':
                pass
            else:
                r = float(data[os.path.splitext(file)[0]][0]/255)
                g = float(data[os.path.splitext(file)[0]][1]/255)
                b = float(data[os.path.splitext(file)[0]][2]/255)
                bpy.ops.import_scene.obj(
                    filepath="{dir}/obj/".format(dir=subjDir) + file)
                mat = bpy.data.materials.new("brainMaterial")
                mat.diffuse_color = (r, g, b, 1)
                bpy.data.objects[os.path.splitext(
                    file)[0]].active_material = mat
                mat.use_nodes = True
                if file[2] == '.':
                    bpy.data.objects[os.path.splitext(
                        file)[0]].parent = bpy.data.objects['Gyri']
                elif file == 'Left-Cerebral-White-Matter.obj' or file == 'Right-Cerebral-White-Matter.obj':
                    bpy.data.objects[os.path.splitext(
                        file)[0]].parent = bpy.data.objects['WhiteMatter']
                else:
                    bpy.data.objects[os.path.splitext(
                        file)[0]].parent = bpy.data.objects['Brain']
                # bpy.data.objects['WhiteMatter'].parent = bpy.data.objects['Brain']
                # bpy.data.objects['Gyri'].parent = bpy.data.objects['Brain']
    else:
        bpy.ops.import_scene.obj(
            filepath="{dir}/obj/".format(dir=subjDir) + 'Left-Cerebral-Cortex.obj')
        mat = bpy.data.materials.new("brainMaterial")
        mat.diffuse_color = (
            float(245 / 255), float(245 / 255), float(245 / 255), 1)
        o = bpy.context.selected_objects[0]
        o.active_material = mat

        bpy.ops.import_scene.obj(
            filepath="{dir}/obj/".format(dir=subjDir) + 'Right-Cerebral-Cortex.obj')
        o = bpy.context.selected_objects[0]
        o.active_material = mat

    bpy.ops.export_scene.fbx(
        filepath="{dir}/{patient}".format(dir=subjDir, patient="reconstruction") + ".fbx")
    bpy.ops.export_scene.gltf(
        export_format="GLB",
        filepath="{dir}/{patient}".format(dir=subjDir,
                                          patient="reconstruction"),
        export_texcoords=False,
        export_normals=False)


if __name__ == "__main__":
    initial_run = sys.argv[5].lower() == 'true'

    main(sys.argv[4], initial_run, sys.argv[6])

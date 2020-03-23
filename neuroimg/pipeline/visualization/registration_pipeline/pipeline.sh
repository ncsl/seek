# #!/bin/bash
mkdir $SUBJECTS_DIR/$1/obj
mkdir $SUBJECTS_DIR/$1/rois
./scripts/aseg2srf.sh $1

mris_convert $SUBJECTS_DIR/$1/surf/lh.pial $SUBJECTS_DIR/$1/surf/lh.pial.asc
mris_convert $SUBJECTS_DIR/$1/surf/rh.pial $SUBJECTS_DIR/$1/surf/rh.pial.asc
mv $SUBJECTS_DIR/$1/surf/lh.pial.asc $SUBJECTS_DIR/$1/surf/lh.pial.srf
mv $SUBJECTS_DIR/$1/surf/rh.pial.asc $SUBJECTS_DIR/$1/surf/rh.pial.srf
./scripts/annot2dpv $SUBJECTS_DIR/$1/label/lh.aparc.annot $SUBJECTS_DIR/$1/label/lh.aparc.annot.dpv
./scripts/annot2dpv $SUBJECTS_DIR/$1/label/rh.aparc.annot $SUBJECTS_DIR/$1/label/rh.aparc.annot.dpv

./scripts/splitsrf $SUBJECTS_DIR/$1/surf/lh.pial.srf $SUBJECTS_DIR/$1/label/lh.aparc.annot.dpv $SUBJECTS_DIR/$1/rois/lh.pial_roi
./scripts/splitsrf $SUBJECTS_DIR/$1/surf/rh.pial.srf $SUBJECTS_DIR/$1/label/rh.aparc.annot.dpv $SUBJECTS_DIR/$1/rois/rh.pial_roi

counter=0001
while [ $counter -le 0035 ]
do
	newCount=$(printf "%04g" $counter)
	labelName=$(cat ./scripts/roiNames.json | jq '.["'$newCount'"]')
	labelName="${labelName#\"}"
	labelName="${labelName%\"}"
	./scripts/srf2obj $SUBJECTS_DIR/$1/rois/lh.pial_roi.$newCount.srf > $SUBJECTS_DIR/$1/obj/lh.$labelName.obj
	./scripts/srf2obj $SUBJECTS_DIR/$1/rois/rh.pial_roi.$newCount.srf > $SUBJECTS_DIR/$1/obj/rh.$labelName.obj
	((counter++))
done


# # Segment the brainstem
# segmentBS.sh $1
# # Segment the thalamic nuclei
# segmentThalamicNuclei.sh $1 $SUBJECTS_DIR
# segmentThalamicNuclei.sh  $1 $SUBJECTS_DIR T2.mgz T2 't2'
# # Segment the hippocampus/amygdala
# segmentHA_T2.sh $1 T2.mgz T2 1
# segmentHA.sh $1
# # On the other hand ThalamicNuclei.v10.T1.FSvoxelSpace.mgz lives in the same voxel space as the other FreeSurfer volumes (e.g., orig.mgz, nu.mgz, aseg.mgz), so you can use it directly to produce masks for further analyses, but its resolution is lower (1 mm vs 0.5 mm).

# # patient, electrodeExport, justCortex
# /usr/local/blender/blender --background --python scripts/sceneCreator.py /$1 True False
#!/bin/bash

# List of labels to be converted if no list is specified
LABLIST="2 3 4 5 7 8 10 11 12 13 14 15 16 17 18 24 26 28 30 31 41 42 43 44 46 47 49 50 51 52 53 54 58 60 62 63 72 77 78 79 80 81 82 85 251 252 253 254 255"

# Prepare a random string to save temporary files
if hash md5 2> /dev/null ; then
    RND0=$(head -n 1 /dev/random | md5)
    elif hash md5sum 2> /dev/null ; then
    RND0=$(head -n 1 /dev/random | md5sum)
fi
RNDSTR=${RND0:0:12}

# Define a function for Ctrl+C as soon as the RNDSTR is defined
trap bashtrap INT
bashtrap()
{
    break ; break
    [[ "${s}" != "" ]] && rm -rf ${SUBJECTS_DIR}/${s}/tmp/${RNDSTR} # try to delete temp files
    exit 1
}

s=$1
# Create directories for temp files and results
mkdir -p ${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}
mkdir -p ${SUBJECTS_DIR}/${s}/ascii

# For each label
for lab in ${LABLIST} ; do
	
	# Label string
	lab0=$(printf %03d ${lab})

	# Pre-tessellate
	echo "==> Pre-tessellating: ${s}, ${lab0}"
	${FREESURFER_HOME}/bin/mri_pretess \
	${SUBJECTS_DIR}/${s}/mri/aseg.mgz ${lab} \
	${SUBJECTS_DIR}/${s}/mri/norm.mgz \
	${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0}_filled.mgz
	
	# Tessellate
	echo "==> Tessellating: ${s}, ${lab0}"
	${FREESURFER_HOME}/bin/mri_tessellate \
	${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0}_filled.mgz \
	${lab} ${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0}_notsmooth
	
	# Smooth
	echo "==> Smoothing: ${s}, ${lab0}"
	${FREESURFER_HOME}/bin/mris_smooth -nw \
	${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0}_notsmooth \
	${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0}
	
    # Convert to ASCII
    echo "==> Converting to ASCII: ${s}, ${lab0}"
    ${FREESURFER_HOME}/bin/mris_convert \
           ${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0} \
           ${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0}.asc
    mv ${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}/aseg_${lab0}.asc \
       ${SUBJECTS_DIR}/${s}/ascii/aseg_${lab0}.srf


	labelName=$(cat ./scripts/LUT.json | jq '.["'$lab'"]')
	labelName="${labelName#\"}"
	labelName="${labelName%\"}"

	./scripts/srf2obj $SUBJECTS_DIR/$1/ascii/aseg_${lab0}.srf > $SUBJECTS_DIR/$1/obj/$labelName.obj

done

rm -rf ${SUBJECTS_DIR}/${s}/tmp/${RNDSTR}

echo "==> Done: ${s}"

exit 0

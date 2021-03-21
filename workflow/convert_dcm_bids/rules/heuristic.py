import os
from heudiconv import lgr as logger

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes


def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where

    allowed template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    # heuristic key for T1w data
    t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item:02d}_T1w')

    # heuristic key for CT data
    ct = create_key('sub-{subject}/{session}/ct/sub-{subject}_{session}_run-{item:02d}_CT')

    info = {t1w: [], ct: []}

    for s in seqinfo:
        """
        The namedtuple `s` contains the following fields:

        * total_files_till_now
        * example_dcm_file
        * series_id
        * dcm_dir_name
        * unspecified2
        * unspecified3
        * dim1
        * dim2
        * dim3
        * dim4
        * TR
        * TE
        * protocol_name
        * is_motion_corrected
        * is_derived
        * patient_id
        * study_description
        * referring_physician_name
        * series_description
        * image_type
        """
        print('The sequence is... ', s)
        logger.info(s)
        if 'mri' in s.dcm_dir_name:
            if 'premri' in s.dcm_dir_name:
                info[t1w].append(s.series_id)
        elif 'ct' in s.dcm_dir_name:
            info[ct].append(s.series_id)
    return info

Semi-Automated Localizing Electrodes Process
============================================

For NCSL specific Readme, check out: `link <extra_docs/contact_localization/localizingelectrodes_instructions.pdf>`_
For general-purpose Readme, check out: `link <extra_docs/contact_localization/localizingelectrodes_instructions.pdf>`_

To only localize contacts using fieldtrip toolbox GUI, or img_pipe GUI, then follow these instructions:


#. Install SPM (preferably 12): https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
#. Install FieldTripToolbox: (http://www.fieldtriptoolbox.org/download/)

Running Localization GUI
------------------------

This assumes you have already ran reconstruction on your T1 MRI and have preprocessed the CT image and downloaded the 
necessary files.


#. 
   Use matlab script to get Voxel/MM coords in CT space

    This requires the user to first have preprocessed the CT scans (and optionally the T1 MRI). 

   .. code-block::

       matlab ./pipeline/contact_localization/matlab/run_localization_fieldtrip.m


    This will run an ~10-15 minute process to have users determine how to localize the channels. Note that
    you will need the corresponding implantation map (i.e. PPT, some image drawn up by clinician, or the implantation knowledge).
    Deep channels (i.e. A1, B1, B'1, etc.) are usually in the brain, while the last channels of
    an electrode are near the skull. 

#. 
   Apply coregistration transform matrix to coords to map to your MRI space.

    TBD

#. 
   (optional) Apply additional affine transformations to make sure your now T1 coordinates
   are in an appropriate coordinate system and/or space, so that when using coords with atlas labels,
   surface files and other T1-extracted image volumes, the coords are in the same 
   language.

    TBD

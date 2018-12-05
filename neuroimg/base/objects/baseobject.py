import os
import warnings


class BaseMeta(object):
    """
    Base class for any metadata object related to the clinical data
    """

    def set_filename(self):
        raise NotImplementedError("Need to implement a function for setting the"
                                  "new filename for the resulting dataset.")

    def create_metadata(self):
        raise NotImplementedError("Need to implement a create_metadata function for creating"
                                  "the metadata for the resulting dataset.")

    @staticmethod
    def format_metadata(metadata):
        if 'aslp' in metadata['filename'].lower():
            metadata['type'] = 'ii asleep'
        if 'aw' in metadata['filename'].lower():
            metadata['type'] = 'ii awake'

        metadata['filename'] = BaseMeta.map_filename(metadata['filename'].lower())

        if 'ii' in metadata['filename']:
            metadata['type'] = 'ii'
        if 'sz' in metadata['filename']:
            metadata['type'] = 'sz'

        return metadata

    @staticmethod
    def map_filename(filename):
        ictal_marks = ['ictal', 'seiz', 'sz']
        ii_marks = ['aslp', 'aw', 'interictal', 'ii', 'inter']

        filename = filename.lower()

        # consolidate markers into one type
        for mark in ictal_marks:
            if mark in filename:
                filename = filename.replace(mark, 'sz')
        for mark in ii_marks:
            if mark in filename:
                filename = filename.replace(mark, 'ii')

        # replace into a fif file format
        filename = filename.replace('.edf', '_raw.fif')
        filename = filename.replace('_0001', '')
        # get rid of spaces
        filename = filename.replace(' ', '')
        # always add _ after type of dataset identifier
        filename = filename.replace('sz', '_sz_')
        filename = filename.replace('ii', '_ii_')

        filename = filename.replace('__', '_')

        # get rid of 'update' in jhu pats
        filename = filename.replace('updated', '')

        return filename

    @staticmethod
    def check_metadata(metadata):
        keys = metadata.keys()
        assert 'filename' in keys
        assert 'patient_id' in keys
        assert 'clinical_center' in keys
        assert all(['onset' in keys, 'termination' in keys])
        assert all(['bad_channels' in keys, 'non_eeg_channels' in keys])

        # check filename
        filename = metadata['filename']
        if os.path.basename(filename) != filename:
            warnings.warn("Filename is not just filename, but a filepath. Make sure to fix\
                for {}!".format(filename))

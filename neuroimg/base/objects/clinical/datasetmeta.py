import os
from edp.base.objects.baseobject import BaseMeta

class DatasetMeta(BaseMeta):
    def __init__(self, dataset_id):
        self.dataset_id = dataset_id

        self.filename = None
        self.equipment = ''
        self.date_of_recording = ''
        self.length_of_recording = -1
        self.number_chans = -1
        self.type = None
        self.note = ""
        self.patient_id = None
        self.onset = None
        self.termination = None
        self.bad_channels = []
        self.non_eeg_channels = []

        '''
        TO BE DEPRECATED
        '''
        self.clinical_center = None
        self.outcome = None
        self.engel_score = None
        self.ez_contacts = []
        self.resect_contacts = []

    def _set_filename(self, filename):
        # ensure that this is just a filename
        filename = os.path.basename(filename).lower()
        return filename

    def set_bad_channels(self, bad_channels):
        bad_channels = [ch.lower() for ch in bad_channels]
        self.bad_channels = bad_channels

    def set_non_eeg_channels(self, non_eeg_channels):
        non_eeg_channels = [ch.lower() for ch in non_eeg_channels]
        self.non_eeg_channels = non_eeg_channels

    def set_dataset_data(self, filename='', equipment='', date_of_recording='', length_of_recording='',
                         number_chans='', type='', note='', onset=[], termination=[], events=[]):
        self.filename = self._set_filename(filename)
        self.equipment = equipment
        self.date_of_recording = date_of_recording
        self.length_of_recording = length_of_recording
        self.number_chans = number_chans
        self.type = type
        self.note = note
        self.onset = onset
        self.termination = termination
        self.events = events

    def set_clinical_data(self, patient_id='', clinical_center='', outcome='', engel_score=None,
                          ez_contacts=[], resect_contacts=[]):
        self.patient_id = patient_id
        self.clinical_center = clinical_center
        self.outcome = outcome
        self.engel_score = engel_score
        self.ez_contacts = ez_contacts
        self.resect_contacts = resect_contacts

    def create_metadata(self):
        jsonobj = {
            'filename': self.filename,
            'equipment': self.equipment,
            'date_of_recording': self.date_of_recording,
            'length_of_recording': self.length_of_recording,
            'number_chans': self.number_chans,
            'type': self.type,
            'note': self.note,
            'patient_id': self.patient_id,
            'dataset_id': self.dataset_id,
            'onset': self.onset,
            'termination': self.termination,
            'bad_channels': self.bad_channels,
            'non_eeg_channels': self.non_eeg_channels,

            'events': self.events,

            # these will be DEPRECATED IN THE FUTURE
            'clinical_center': self.clinical_center,
            'outcome': self.outcome,
            'engel_score': self.engel_score,
            'ez_elecs': self.ez_contacts,
            'resect_elecs': self.resect_contacts
        }
        super(DatasetMeta, self).check_metadata(jsonobj)

        jsonobj = super(DatasetMeta, self).format_metadata(jsonobj)
        return jsonobj

import mne
import numpy as np

from edp.base.objects.dataset.basedataobject import BaseDataset
from edp.format.base.elecs import Contacts
from edp.utils.scalprecording import create_mne_topostruct_from_numpy


class TimeSeries(BaseDataset):
    """
    The base class for our time series data.

    All time series are assumed to be in [C x T] shape and use the Contacts
    data structure to handle all contact level functionality.

    TimeSeries -> Metadata:
                    - Contacts
                    - metadata json object

    Properties:
    - contact_tuple_list
    - electrodes
    - chanlabels
    - xyz_coords
    - ncontacts
    - len_secs
    """

    def __init__(self, mat, metadata):
        times = np.arange(mat.shape[1]).astype(int)
        super(TimeSeries, self).__init__(mat=mat, times=times)

        self.metadata = metadata
        if self.mat.ndim > 2:
            raise ValueError("Time series can not have > 2 dimensions right now."
                             "We assume [C x T] shape, channels x time. ")

        # extract metadata for this time series
        self.extract_metadata()

    def get_metadata(self):
        return self.metadata

    @property
    def length_of_recording(self):
        return len(self.times)

    @property
    def len_secs(self):
        """
        Determines the length of the recording in seconds
        :return:
        """
        return self.length_of_recording / float(self.samplerate)

    def extract_metadata(self):
        self.samplerate = self.metadata['samplerate']
        self.patient_id = self.metadata['patient_id']
        self.dataset_id = self.metadata['dataset_id']
        self.contacts_list = self.metadata['chanlabels']
        # self.clinonsetlabels = self.metadata['ez_hypo_contacts']

        try:
            self.cezlobe = self.metadata['cezlobe']
        except:
            self.cezlobe = []

        try:
            self.clinonsetlabels = self.metadata['clinezelecs']
            self.clinonsetlabels = self.metadata['ez_hypo_contacts']
        except Exception as e:
            print(e)
            self.clinonsetlabels = []

        try:
            self.outcome = self.metadata['outcome']
        except Exception as e:
            print("Error in extracting metadata: ", e)
            self.outcome = ''

        try:
            self.montage = self.metadata['montage']
        except:
            self.montage = ''

        try:
            self.clinicaldifficulty = self.metadata['clinical_difficulty']
        except:
            self.clinicaldifficulty = ''

        try:
            self.clinicalmatching = self.metadata['clinical_match']
        except:
            self.clinicalmatching = ''

        self.onsetind = self.metadata['onsetind']
        self.offsetind = self.metadata['offsetind']

        # convert channel labels into a Contacts data struct
        self.contacts = Contacts(self.contacts_list, require_matching=False)

    def set_bipolar(self, chanlabels=[]):
        """
        Set bipolar montage for the eeg data time series.

        :param chanlabels:
        :return:
        """
        # extract bipolar reference scheme from contacts data structure
        self._bipolar_inds = self.contacts.set_bipolar(chanlabels=chanlabels)

        # set the time series to be bipolar
        self.mat = self.mat[self._bipolar_inds[:, 1], :] - \
                   self.mat[self._bipolar_inds[:, 0], :]
        self.metadata['chanlabels'] = self.chanlabels

    def set_common_avg_ref(self):
        """
        Set a common average referencing scheme.

        :return:
        """
        self.ref_signal = np.mean(self.mat, axis=0)
        self.mat = self.mat - self.ref_signal

    def set_reference_signal(self, ref_signal):
        """
        Set a custom reference signal to reference all signals by.

        :param ref_signal: (np.ndarray) [Tx1] vector that is subtracted from all signals
        in our dataset.
        :return:
        """
        self.ref_signal = ref_signal
        self.mat = self.mat - self.ref_signal

    def filter_data(self, linefreq, samplerate):
        """
        Filters the time series data according to the line frequency (notch) and
        sampling rate (band pass filter).

        :param linefreq:
        :param samplerate:
        :return:
        """
        # the bandpass range to pass initial filtering
        freqrange = [0.5]
        freqrange.append(samplerate // 2 - 1)

        # the notch filter to apply at line freqs
        linefreq = int(linefreq)  # LINE NOISE OF HZ
        assert linefreq == 50 or linefreq == 60

        # initialize the line freq and its harmonics
        freqs = np.arange(linefreq, samplerate // 2, linefreq)
        freqs = np.delete(freqs, np.where(freqs > samplerate // 2)[0])

        # run bandpass filter and notch filter
        self.mat = mne.filter.filter_data(self.mat,
                                          sfreq=samplerate,
                                          l_freq=freqrange[0],
                                          h_freq=freqrange[1],
                                          # pad='reflect',
                                          verbose=False
                                          )
        self.mat = mne.filter.notch_filter(self.mat,
                                           Fs=samplerate,
                                           freqs=freqs,
                                           verbose=False
                                           )

    def mask_indices(self, mask_inds):
        self.mat = self.mat[mask_inds, :]
        self.contacts.mask_contact_indices(mask_inds)

    def create_raw_struct(self, remove=True):
        eegts = self.get_data()
        chanlabels = self.chanlabels
        metadata = self.get_metadata()
        samplerate = self.samplerate

        # if montage not in mne.
        # create the info struct
        info = mne.create_info(ch_names=chanlabels.tolist(), ch_types=['eeg'] * len(chanlabels), sfreq=samplerate)
        # create the raw object
        mne_rawstruct = mne.io.RawArray(data=eegts, info=info)
        return mne_rawstruct

    def trim_aroundonset(self, offset_sec=20, mat=None):
        """
        Trims dataset to have (seconds) before/after onset/offset.

        If there is no offset, then just takes it offset after onset.

        :param offset:
        :return:
        """
        # get 30 seconds before/after
        offsetind = int(offset_sec * self.samplerate)
        preindex = self.onsetind - offsetind
        postindex = self.onsetind + offsetind
        interval = (preindex, postindex)

        if mat is None:
            mat, times = self.trim_dataset(interval=interval)
        else:
            mat = mat[..., preindex:postindex]
            times = self.times[preindex:postindex]

        self.mat = mat
        self.times = times

        return mat, offsetind

class EEGTs(TimeSeries):
    def __init__(self, mat, metadata, montage):
        super(EEGTs, self).__init__(mat, metadata)
        self.montage = montage

    def __str__(self):
        return "{} {} EEG mat ({}) " \
               "{} seconds".format(self.patient_id, self.dataset_id, self.mat.shape, self.len_secs)

    def __repr__(self):
        return "{} {} EEG mat ({}) " \
               "{} seconds".format(self.patient_id, self.dataset_id, self.mat.shape, self.len_secs)

    def apply_montage(self, montage=None):
        montage = self.montage
        if self.montage is None:
            montage = 'standard_1020'

        # create a mne raw data structure with the montage applied
        mne_rawstruct = create_mne_topostruct_from_numpy(self.mat,
                                                         self.chanlabels,
                                                         self.samplerate,
                                                         montage=montage)
        return mne_rawstruct

    def set_bipolar(self, chanlabels=[]):
        raise RuntimeError("It is not recommended to set bipolar contacts"
                           "on scalp EEG! Should use for SEEG and ECoG.")

    def create_raw_struct(self, remove=True):
        eegts = self.get_data()
        chanlabels = self.chanlabels
        metadata = self.get_metadata()
        samplerate = self.samplerate

        # gets the best montage
        best_montage = self.get_best_matching_montage(chanlabels)

        # gets channel indices to keep
        montage_chs = chanlabels
        montage_data = eegts
        if remove:
            montage_inds = self.get_montage_channel_indices(best_montage, chanlabels)
            montage_chs = chanlabels[montage_inds]
            other_inds = [idx for idx, ch in enumerate(chanlabels) if ch not in montage_chs]
            montage_data = eegts[montage_inds, :]

            print("Removed these channels: ", chanlabels[other_inds])

        mne_rawstruct = create_mne_topostruct_from_numpy(montage_data, montage_chs, samplerate, montage=best_montage)
        return mne_rawstruct

    def get_montage_channel_indices(self, montage_name, chanlabels):
        # read in montage and strip channel labels not in montage
        montage = mne.channels.read_montage(montage_name)
        montage_chs = [ch.lower() for ch in montage.ch_names]

        # get indices to keep
        to_keep_inds = [idx for idx, ch in enumerate(chanlabels) if ch in montage_chs]

        return to_keep_inds

    def get_best_matching_montage(self, chanlabels):
        """
        Get the best matching montage with respect to this montage
        :param chanlabels:
        :return:
        """
        montages = mne.channels.get_builtin_montages()
        best_montage = None
        best_montage_score = 0

        for montage_name in montages:
            # read in standardized montage
            montage = mne.channels.read_montage(montage_name)

            # get the channels and score for this montage wrt channels
            montage_score = 0
            montage_chs = [ch.lower() for ch in montage.ch_names]

            # score this montage
            montage_score = len([ch for ch in chanlabels if ch in montage_chs])

            if montage_score > best_montage_score:
                best_montage = montage_name
                best_montage_score = montage_score

        return best_montage


class SEEGTs(TimeSeries):
    def __init__(self, mat, metadata):
        super(SEEGTs, self).__init__(mat, metadata)

    def __str__(self):
        return "{} {} SEEG mat ({}) " \
               "{} seconds".format(self.patient_id, self.dataset_id, self.mat.shape, self.len_secs)

    def get_hemisphere_contacts(self, hemisphere):
        if hemisphere == 'r':
            contacts = [ch for idx, ch in enumerate(self.chanlabels) if "'" not in ch]
        elif hemisphere == 'l':
            contacts = [ch for idx, ch in enumerate(self.chanlabels) if "'" in ch]
        return contacts

class ECOGTs(TimeSeries):
    def __init__(self):
        super(ECOGTs, self).__init__()

    def __str__(self):
        return "{} {} ECoG mat ({}) " \
               "{} seconds".format(self.patient_id, self.dataset_id, self.mat.shape, self.len_secs)

    def get_grid_contacts(self):
        grid_contacts = [ch for idx, ch in enumerate(self.chanlabels) if ch.startswith('g')]

        return grid_contacts
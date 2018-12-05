import re

import mne
import numpy as np

from edp.base.config.dataconfig import Freqbands
from edp.base.objects.dataset.basedataobject import BaseDataset
from edp.format.base.elecs import Contacts
from edp.utils.post_process import PostProcess
from edp.utils.utils import compute_timepoints, map_to_win


class Result(BaseDataset):
    def __init__(self, mat, metadata):
        self.metadata = metadata

        # extract metadata for this time series
        self.extract_metadata()

        super(Result, self).__init__(mat=mat, times=self.timepoints[:, 1])

    def get_metadata(self):
        return self.metadata

    def compute_montage_groups(self):
        from mne.selection import _divide_to_regions
        rawinfo = mne.create_info(ch_names=self.chanlabels.tolist(),
                                  ch_types='eeg',
                                  sfreq=self.samplerate,
                                  montage=self.montage)

        # get channel groups - hashmap of channel indices
        ch_groups = _divide_to_regions(rawinfo, add_stim=False)

        # add this to class object but with lower-cased keys
        self.ch_groups = {}
        for k, v in ch_groups.items():
            self.ch_groups[k.lower()] = v

        # get the indices of the cez lobe
        self.cezlobeinds = []
        for lobe in self.cezlobe:
            self.cezlobeinds.extend(self.ch_groups[lobe])
        self.oezlobeinds = [ind for ind in range(len(self.chanlabels)) if ind not in self.cezlobeinds]

        print(len(self.chanlabels))

    def extract_metadata(self):
        self.samplerate = self.metadata['samplerate']
        self.patient_id = self.metadata['patient_id']
        self.dataset_id = self.metadata['dataset_id']
        self.contacts_list = self.metadata['chanlabels']
        self.samplepoints = np.array(self.metadata['samplepoints'])

        try:
            self.cezlobe = self.metadata['cezlobe']
        except:
            self.cezlobe = []

        # self.modality = self.metadata['modality']
        try:
            self.outcome = self.metadata['outcome']
            self.clinonsetlabels = self.metadata['ez_hypo_contacts']
            self.semiology = self.metadata['seizure_semiology']
        except:
            self.outcome = ''
            self.clinonsetlabels = []
            self.semiology = []

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

        try:
            # comment out
            self.modality = self.metadata['modality']
        except Exception as e:
            self.metadata['modality'] = 'ieeg'
            self.modality = self.metadata['modality']
            print("Error in extracting metadata: ", e)

        self.onsetwin = self.metadata['onsetwin']
        self.offsetwin = self.metadata['offsetwin']
        self.winsize = self.metadata['winsize']
        self.stepsize = self.metadata['stepsize']
        self.reference = self.metadata['reference']

        if 'radius' in self.metadata.keys():
            self.radius = self.metadata['radius']
        if 'perturbtype' in self.metadata.keys():
            self.perturbtype = self.metadata['perturbtype']

        # compute time points
        self.timepoints = self.samplepoints.astype(int) / self.samplerate

        # convert channel labels into a Contacts data struct
        if self.reference == 'bipolar' or self.modality == 'scalp':
            self.contacts = Contacts(self.contacts_list, require_matching=False)
        else:
            self.contacts = Contacts(self.contacts_list)

    def mask_channels(self):
        badchannels = self.metadata['bad_channels']
        noneegchannels = self.metadata['non_eeg_channels']

        maskinds = []
        for chlist in [badchannels, noneegchannels]:
            removeinds = self.remove_channels(chlist)
            maskinds.extend(removeinds)

        maskinds = list(set(maskinds))
        nonmaskinds = [ind for ind in range(len(self.chanlabels)) if ind not in maskinds]

        # apply to relevant data
        self.mat = self.mat[nonmaskinds, ...]
        self.contacts.mask_contact_indices(nonmaskinds)

    def trim_aroundonset(self, offset_sec=20, mat=None):
        """
        Trims dataset to have (seconds) before/after onset/offset.

        If there is no offset, then just takes it offset after onset.

        :param offset:
        :return:
        """
        # get 30 seconds before/after
        offsetind = int(offset_sec * self.samplerate)
        offsetwin = map_to_win(offsetind, samplepoints=self.samplepoints)
        preindex = self.onsetwin - offsetwin
        postindex = self.onsetwin + offsetwin
        # if preindex.size > 1:
        #     preindex = preindex[0]
        # if postindex.size > 1:
        #     postindex = postindex[0]

        interval = (preindex, postindex)
        print(interval)
        if mat is None:
            mat, times = self.trim_dataset(interval=interval)
        else:
            mat = mat[..., preindex:postindex]
            times = self.times[preindex:postindex]

        self.mat = mat
        self.times = times
        print(self.times.shape)
        print(self.mat.shape)
        return mat, offsetwin

    def trim_aroundseizure(self, offset_sec=20, mat=None, apply_postprocess=True):
        """
        Trims dataset to have (seconds) before/after onset/offset.

        If there is no offset, then just takes it offset after onset.

        :param offset_sec:
        :return:
        """
        # parameters of the resampling length
        ICTAL_LEN = 100
        PREICTAL_LEN = POSTICTAL_LEN = 30

        # get 30 seconds before/after
        offsetind = int(offset_sec * self.samplerate)
        offsetwin = map_to_win(offsetind, samplepoints=self.samplepoints)

        # set the first index cutoff
        if self.onsetwin is not None:
            preindex = self.onsetwin - offsetwin
            if preindex < 0:
                preindex = 0
        else:
            preindex = 0

        # set the postindex
        if self.offsetwin is not None and not np.isnan(self.offsetwin):
            postindex = self.offsetwin + offsetwin
        else:
            postindex = self.onsetwin + offsetwin * 2

        interval = (preindex, postindex)
        # print("Trimming dataset:", mat.shape, preindex, postindex, offsetwin)
        if mat is None:
            mat, times = self.trim_dataset(interval=interval)
        else:
            mat = mat[..., preindex:postindex]
            times = self.times[preindex:postindex]

        print(self.onsetwin, self.offsetwin)
        print(mat.shape, offsetwin)

        if apply_postprocess:
            # warp the ictal region
            warped_ictal = PostProcess.resample_mat(mat[..., offsetwin:-offsetwin],
                                                    ICTAL_LEN)

            # get the preictal region
            preictal_mat = PostProcess.resample_mat(mat[..., :offsetwin],
                                                    PREICTAL_LEN)

            # get the postictal region
            postictal_mat = PostProcess.resample_mat(mat[..., offsetwin:],
                                                     POSTICTAL_LEN)

            # concatenate preictal, ictal and postictal warped region
            mat = np.concatenate((preictal_mat, warped_ictal, postictal_mat), axis=1)

            # create new times array that doens't matter what original scale of time was.
            times = np.arange(mat.shape[1])

        self.mat = mat
        self.times = times

        return mat

    def expand_bipolar_chans(self, ch_list):
        if ch_list == []:
            return None

        ch_list = [a.replace("’", "'") for a in ch_list]
        new_list = []
        for string in ch_list:
            if not string.strip():
                continue

            # A1-10
            match = re.match("^([a-z]+)([0-9]+)-([0-9]+)$", string)
            if match:
                name, fst_idx, last_idx = match.groups()

                new_list.extend([name + str(fst_idx), name + str(last_idx)])

        return new_list

    def make_onset_labels_bipolar(self, clinonsetlabels):
        added_ch_names = []
        for ch in clinonsetlabels:
            # A1-10
            match = re.match("^([a-z]+)([0-9]+)$", ch)
            if match:
                name, fst_idx = match.groups()
            added_ch_names.append(name + str(int(fst_idx) + 1))

        clinonsetlabels.extend(added_ch_names)
        clinonsetlabels = list(set(clinonsetlabels))
        return clinonsetlabels


class LtvResult(Result):
    """
    A class wrapper for linear-time varying model result. We assume a structure of windows over time with a model
    imposed on the data as a function of that sliding window. This also allows overlapping windows and as a consequence
    overlapping models.

    """

    def __init__(self, adjmats, metadata):
        if adjmats.ndim != 3:
            raise AttributeError("The passed in ltv model needs to have 3 dimensions!")

        if adjmats.shape[1] == adjmats.shape[2]:
            if adjmats.shape[0] == adjmats.shape[1]:
                raise RuntimeError("Need to implement how to roll back axis for LTV result here!")

            # make sure the time axis is rolled to the back
            adjmats = np.moveaxis(adjmats, 0, -1)

        super(LtvResult, self).__init__(adjmats, metadata)

        self.json_fields = [
            'onsetwin',
            'offsetwin',
            'resultfilename',
            'winsize',
            'stepsize',
        ]

    def __str__(self):
        return "{} {} LTV Model {}".format(self.patient_id,
                                           self.dataset_id,
                                           self.shape)


class PerturbationResult(Result):

    def __init__(self, pertmat, metadata):
        super(PerturbationResult, self).__init__(pertmat, metadata)

        if pertmat.ndim != 2:
            raise AttributeError("The passed in perturbation model needs to have only 2 dimensions!")

        self.json_fields = [
            'onsetwin',
            'offsetwin',
            'resultfilename',
            'winsize',
            'stepsize',
            'radius',
            'perturbtype',
        ]

    def __str__(self):
        return "{} {} Min-Norm Perturbation Model {}".format(self.patient_id,
                                                             self.dataset_id,
                                                             self.shape)


class FreqModelResult(Result):

    def __init__(self, mat, metadata, phasemat=[], freqbands=Freqbands):
        super(FreqModelResult, self).__init__(mat, metadata)
        self.phasemat = phasemat

        # initialize the frequency bands
        self.freqbands = freqbands
        # set the frequency index
        self.freqindex = 1

        # create a copy of the matrix
        self.buffmat = self.mat.copy()

    def get_phase(self):
        return self.phasemat

    def get_freqs(self):
        return self.freqs

    def extract_metadata(self):
        self.samplerate = self.metadata['samplerate']
        self.patient_id = self.metadata['patient_id']
        self.dataset_id = self.metadata['dataset_id']
        self.contacts_list = self.metadata['chanlabels']
        self.samplepoints = np.array(self.metadata['samplepoints'])

        self.onsetwin = self.metadata['onsetwin']
        self.offsetwin = self.metadata['offsetwin']
        self.winsize = self.metadata['winsize']
        self.stepsize = self.metadata['stepsize']

        self.freqs = self.metadata['freqs']

        self.modality = self.metadata['modality']
        # self.outcome = self.metadata['outcome']
        # self.clinonsetlabels = self.metadata['clinezelecs']
        try:
            self.clinonsetlabels = self.metadata['clinezelecs']
        except Exception as e:
            print("Error in extracting metadata: ", e)
            self.clinonsetlabels = []

        try:
            self.outcome = self.metadata['outcome']
        except:
            self.outcome = ''

        try:
            self.montage = self.metadata['montage']
        except:
            self.montage = ''

        try:
            self.cezlobe = self.metadata['cezlobe']
        except:
            self.cezlobe = []

        try:
            # comment out
            self.modality = self.metadata['modality']
        except Exception as e:
            self.metadata['modality'] = 'ieeg'
            self.modality = self.metadata['modality']
            print("Error in extracting metadata: ", e)

        # compute time points
        self.timepoints = compute_timepoints(self.samplepoints.ravel()[-1],
                                             self.winsize,
                                             self.stepsize,
                                             self.samplerate)

        # convert channel labels into a Contacts data struct -> don't require matching between contacts
        self.contacts = Contacts(self.contacts_list, require_matching=False)

    def __str__(self):
        return "{} {} Frequency Power Model {}".format(self.patient_id,
                                                       self.dataset_id,
                                                       self.shape)

    def _computefreqindices(self, freqs, freqband):
        """
        Compute the frequency indices for this frequency band [lower, upper].

        freqs = list of frequencies
        freqband = [lowerbound, upperbound] frequencies of the
                frequency band
        """
        lowerband = freqband[0]
        upperband = freqband[1]

        # get indices where the freq bands are put in
        freqbandindices = np.where(
            (freqs >= lowerband) & (freqs < upperband))
        freqbandindices = [freqbandindices[0][0], freqbandindices[0][-1]]
        return freqbandindices

    def compress_freqbands(self, freqbands=Freqbands):
        # ensure power is absolute valued
        power = np.abs(self.mat)

        # create empty binned power
        power_binned = np.zeros(shape=(power.shape[0],
                                       len(freqbands),
                                       power.shape[2]))

        for idx, freqband in enumerate(freqbands):
            # compute the freq indices for each band
            freqbandindices = self._computefreqindices(self.freqs, freqband.value)

            # Create an empty array = C x T (frequency axis is compresssed into 1 band)
            # average between these two indices
            power_binned[:, idx, :] = np.mean(power[:, freqbandindices[0]:freqbandindices[1] + 1, :],
                                              axis=1)
        self.mat = power_binned
        self.freqbands = freqbands
        self.freqbandslist = list(freqbands)

    def format_dataset(self, result, freqindex=None, apply_trim=True, is_scalp=False):
        # number of seconds to offset trim dataset by
        offset_sec = 10
        # threshold level
        threshlevel = 0.7

        # get channels separated data by cez and others
        if is_scalp:
            # compute montage group
            result.compute_montage_groups()

            print(result.get_data().shape)
            # print(result.chanlabels)
            # print(result.cezlobe)
            # print(len(result.cezlobeinds), len(result.oezlobeinds))

            # partition into cir and oir
            clinonset_map = result.get_data()[result.cezlobeinds,...]
            others_map = result.get_data()[result.oezlobeinds,...]
        else:
            clinonset_map = result.get_data()[result.cezinds]
            others_map = result.get_data()[result.oezinds]

        if freqindex is not None:
            clinonset_map = clinonset_map[:,freqindex,:]
            others_map = others_map[:, freqindex, :]
            # print(clinonset_map.shape, others_map.shape)

        if apply_trim:
            # trim dataset in time
            clinonset_map = result.trim_aroundseizure(offset_sec=offset_sec, mat=clinonset_map)
            others_map = result.trim_aroundseizure(offset_sec=offset_sec, mat=others_map)

        clinonset_map = np.mean(clinonset_map, axis=0)
        others_map = np.mean(others_map, axis=0)

        return clinonset_map, others_map


class FragilityModelResult(Result):

    def __init__(self, ltvmodel, pertmodel, metadata):
        self.ltvmodel = ltvmodel
        self.pertmodel = pertmodel
        self.mat = FragilityModelResult.compute_fragilitymetric(pertmodel.mat)
        self.minmax_fragmat = FragilityModelResult.compute_minmaxfragilitymetric(pertmodel.mat)

        super(FragilityModelResult, self).__init__(self.mat, metadata)

        self.buffmat = self.mat.copy()

    def __str__(self):
        return "{} {} Fragility Model {}".format(self.patient_id,
                                                 self.dataset_id,
                                                 self.shape)

    def __repr__(self):
        return str(self)

    @staticmethod
    def compute_fragilitymetric(minnormpertmat):
        # get dimensions of the pert matrix
        N, T = minnormpertmat.shape
        # assert N < T
        fragilitymat = np.zeros((N, T))
        for icol in range(T):
            fragilitymat[:, icol] = (np.max(minnormpertmat[:, icol]) - minnormpertmat[:, icol]) / \
                                    np.max(minnormpertmat[:, icol])
        return fragilitymat

    @staticmethod
    def compute_minmaxfragilitymetric(minnormpertmat):
        # get dimensions of the pert matrix
        N, T = minnormpertmat.shape
        # assert N < T
        minmax_fragilitymat = np.zeros((N, T))

        # get the min/max for each column in matrix
        minacrosstime = np.min(minnormpertmat, axis=0)
        maxacrosstime = np.max(minnormpertmat, axis=0)

        # normalized data with minmax scaling
        minmax_fragilitymat = -1 * np.true_divide((minnormpertmat - np.matlib.repmat(maxacrosstime, N, 1)),
                                                  np.matlib.repmat(maxacrosstime - minacrosstime, N, 1))
        return minmax_fragilitymat

    def apply_thresholding_smoothing(self, threshold, mat=None):
        if mat is None:
            mat = self.mat.copy()
            mat[mat < threshold] = 0
            self.mat = mat
        else:
            mat = mat
            mat[mat < threshold] = 0
        return mat

    def apply_time_warping(self):
        pass

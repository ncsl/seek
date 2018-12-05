import numpy as np
from natsort import order_by_index

import mne
from edp.utils.scalprecording import create_mne_topostruct_from_numpy


class BaseDataset(object):
    def __init__(self, mat, times):
        self.mat = mat
        self.times = times

        self.bufftimes = self.times.copy()
        self.buffmat = self.mat.copy()

    def __len__(self):
        return self.mat.shape[1]

    def load_contacts_regs(self, contact_regs, atlas=''):
        self.contacts.load_contacts_regions(contact_regs)
        self.atlas = atlas

    def load_chanxyz(self, chanxyz, coordsystem="T1MRI"):
        """
         Load in the channel's xyz coordinates.

        :param chanxyz:
        :param coordsystem:
        :return:
        """
        if len(chanxyz) != self.ncontacts:
            raise RuntimeError("In loading channels xyz, chanxyz needs to be"
                               "of equal length as the number of contacts in dataset! "
                               "There is a mismatch chanxyz={} vs "
                               "dataset.ncontacts={}".format(
                len(chanxyz), self.ncontacts
            ))
        self.contacts.load_contacts_xyz(chanxyz)
        self.coordsystem = coordsystem

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

    def create_raw_struct(self, usebestmontage=False):
        eegts = self.get_data()
        chanlabels = self.chanlabels
        samplerate = self.samplerate

        # gets the best montage
        best_montage = self.get_best_matching_montage(chanlabels)

        # gets channel indices to keep
        montage_chs = chanlabels
        montage_data = eegts
        if usebestmontage:
            montage_inds = self.get_montage_channel_indices(best_montage, chanlabels)
            montage_chs = chanlabels[montage_inds]
            other_inds = [idx for idx, ch in enumerate(chanlabels) if ch not in montage_chs]
            montage_data = eegts[montage_inds, :]

            print("Removed these channels: ", chanlabels[other_inds])

        mne_rawstruct = create_mne_topostruct_from_numpy(montage_data, montage_chs,
                                                         samplerate, montage=best_montage)
        return mne_rawstruct

    def get_montage_channel_indices(self, montage_name, chanlabels):
        # read in montage and strip channel labels not in montage
        montage = mne.channels.read_montage(montage_name)
        montage_chs = [ch.lower() for ch in montage.ch_names]

        # get indices to keep
        to_keep_inds = [idx for idx, ch in enumerate(chanlabels) if ch in montage_chs]

        return to_keep_inds

    def reset(self):
        self.mat = self.buffmat.copy()
        self.times = self.bufftimes.copy()
        # self.chanlabels = self.buffchanlabels.copy()

    @property
    def shape(self):
        return self.mat.shape

    @property
    def contact_tuple_list(self):
        return [str(a) + str(b) for a, b in self.contacts_list]

    @property
    def electrodes(self):
        # use a dictionary to store all electrodes
        return self.contacts.electrodes

    @property
    def chanlabels(self):
        return np.array(self.contacts.chanlabels)

    @property
    def xyz_coords(self):
        return self.contacts.xyz

    @property
    def ncontacts(self):
        return len(self.chanlabels)

    @property
    def cezinds(self):
        clinonsetinds = [i for i, ch in enumerate(self.chanlabels) if
                         ch in self.clinonsetlabels]
        return clinonsetinds

    @property
    def oezinds(self):
        oezinds = [i for i, ch in enumerate(self.chanlabels) if
                   ch not in self.clinonsetlabels]
        return oezinds

    def natsort_contacts(self):
        """
        Sort out the time series by its channel labels, so they go in
        a natural ordering.

        A1,A2, ..., B1, B2, ..., Z1, Z2, ..., A'1, A'2, ...

        :return:
        """
        print("Trying to sort naturally contacts in result object")

        self.buffchanlabels = self.chanlabels.copy()
        # pass
        natinds = self.contacts.natsort_contacts()
        self.mat = np.array(order_by_index(self.mat, natinds))
        self.mat = np.array(self.mat)
        self.metadata['chanlabels'] = self.chanlabels

    def get_data(self):
        return self.mat

    def get_channel_data(self, name, interval=(None, None)):
        idx = self.chanlabels.index(name)
        tid1, tid2 = self.interval_to_index(interval)
        return self.mat[idx, tid1:tid2]

    def remove_channels(self, toremovechans):
        removeinds = [ind for ind, ch in enumerate(self.chanlabels) if ch in toremovechans]
        return removeinds

    def split_cez_oez(self, dataset, cez_inds, oez_inds):
        return dataset[cez_inds, :], dataset[oez_inds, :]

    def trim_dataset(self, interval=(None, None)):
        """
        Trims dataset to have (seconds) before/after onset/offset.

        If there is no offset, then just takes it offset after onset.

        :param offset:
        :return:
        """
        tid1, tid2 = self.interval_to_index(interval)
        self.mat = self.mat[..., tid1:tid2]
        self.times = self.times[tid1:tid2]

        return self.mat, self.times

    def interval_to_index(self, interval):
        tid1, tid2 = 0, -1
        if interval[0] is not None:
            if interval[0] < self.times[0]:
                tid1 = 0
            else:
                tid1 = np.argmax(self.times >= interval[0])
        if interval[1] is not None:
            if interval[1] > self.times[-1]:
                print(self.times[-1], interval)
                return -1
            else:
                tid2 = np.argmax(self.times >= interval[1])
        return (tid1, tid2)

    def time(self, interval=(None, None)):
        tid1, tid2 = self.interval_to_index(interval)
        return self.times[tid1:tid2]

    def apply_moving_avg_smoothing(self, window_size):
        def movingaverage(interval, window_size):
            window = np.ones(int(window_size)) / float(window_size)
            return np.convolve(interval, window, 'same')

        mat = self.mat.copy()

        # apply moving average filter to smooth out stuff
        smoothed_mat = np.array([movingaverage(x, window_size=window_size) for x in mat])

        # realize that moving average messes up the ends of the window
        smoothed_mat = smoothed_mat[:, window_size // 2:-window_size // 2]

        self.mat = smoothed_mat
        return smoothed_mat

    def apply_gaussian_kernel_smoothing(self, window_size):
        from edp.utils.post_process import PostProcess

        mat = self.mat.copy()

        # apply moving average filter to smooth out stuff
        smoothed_mat = np.array([PostProcess.smooth_kernel(x, window_len=window_size) for x in mat])

        # realize that moving average messes up the ends of the window
        # smoothed_mat = smoothed_mat[:, window_size // 2:-window_size // 2]

        self.mat = smoothed_mat
        return smoothed_mat

    def compute_metric_window(self, winsize, stepsize, metric='avg'):
        """
        Method to compute the average mat value in a specified window.
        :return:
        """
        from edp.utils.utils import compute_samplepoints

        samplepoints = compute_samplepoints(winsize, stepsize, self.mat.shape[1])
        if metric == 'avg':
            pass
        elif metric == 'max':
            pass
        pass

    def format_dataset(self, result, apply_trim=True, apply_thresh=True, is_scalp=False):
        # number of seconds to offset trim dataset by
        offset_sec = 10
        # threshold level
        threshlevel = 0.7

        # get channels separated data by cez and others
        if is_scalp:
            # compute montage group
            result.compute_montage_groups()

            print(result.get_data().shape)
            print(result.chanlabels)
            print(result.cezlobe)
            print(len(result.cezlobeinds), len(result.oezlobeinds))

            # partition into cir and oir
            clinonset_map = result.get_data()[result.cezlobeinds]
            others_map = result.get_data()[result.oezlobeinds]
        else:
            clinonset_map = result.get_data()[result.cezinds]
            others_map = result.get_data()[result.oezinds]

        if apply_trim:
            # trim dataset in time
            clinonset_map = result.trim_aroundseizure(offset_sec=offset_sec, mat=clinonset_map)
            others_map = result.trim_aroundseizure(offset_sec=offset_sec, mat=others_map)

        if apply_thresh:
            # apply thresholding
            clinonset_map = result.apply_thresholding_smoothing(threshold=threshlevel, mat=clinonset_map)
            others_map = result.apply_thresholding_smoothing(threshold=threshlevel, mat=others_map)

        clinonset_map = np.mean(clinonset_map, axis=0)
        others_map = np.mean(others_map, axis=0)

        return clinonset_map, others_map

"""Script to generate <subjID>_setup.mat file for fieldtrip."""
from pathlib import Path
import numpy as np

from mne_bids import BIDSPath, read_raw_bids, get_entity_vals
from scipy.io import savemat


def create_elec_labels_subject(bids_path, output_fpath):
    # read in the actual file
    raw = read_raw_bids(bids_path)

    # get the channels
    ch_names = raw.ch_names

    # save this to a output file
    elec_dict = {
        'elec_name': np.array(ch_names).astype(np.object)
    }
    Path(output_fpath).parent.mkdir(exist_ok=True, parents=True)
    savemat(output_fpath, elec_dict)


if __name__ == '__main__':
    root = Path('/Users/adam2392/Dropbox/epilepsy_bids/')
    # get the runs for this subject
    subjects = get_entity_vals(root, 'subject')
    session = 'presurgery'
    extension = '.vhdr'

    output_path = root / 'sourcedata' / 'electrodes localized' / 'setup'

    for subject in subjects:
        if not any([x in subject for x in ['la', 'nl', 'tvb']]):
            continue
        if subject in ['la00']:
            continue

        # create the BIDS path
        bids_path = BIDSPath(subject=subject, session=session,
                             root=root, extension=extension)

        # get matches
        data_paths = bids_path.match()
        output_fpath = output_path / f'sub-{subject}_setup.mat'
        print(data_paths)
        create_elec_labels_subject(bids_path=data_paths[0], output_fpath=output_fpath)
        # break
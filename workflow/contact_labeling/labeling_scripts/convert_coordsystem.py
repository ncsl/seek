import argparse
import logging

from seek_localize import read_dig_bids, convert_coord_space, write_dig_bids, convert_coord_units

logger = logging.getLogger(__name__)


def main(root, elecs_fname, output_elecs_fname, to_coord_frame, subjects_dir):
    sensors = read_dig_bids(elecs_fname, root=root)

    # convert to voxels
    sensors_vox = convert_coord_units(sensors, to_unit='voxel')

    # convert coordinate frame
    sensors_new = convert_coord_space(sensors_vox, to_frame=to_coord_frame,
                                      subjects_dir=subjects_dir)
    ch_names = sensors_new.ch_names

    # convert coordinate space to xyz
    if sensors_new.coord_unit == 'voxel':
        sensors_new = convert_coord_units(sensors_new, to_unit='mm')

    # save to output specified location on disc
    logger.info(f'saving to... {output_elecs_fname}')
    write_dig_bids(output_elecs_fname, root=root,
                   ch_names=ch_names,
                   ch_coords=sensors_new.get_coords(),
                   coord_system=sensors_new.coord_system,
                   unit=sensors_new.coord_unit,
                   intended_for=sensors_new.intended_for)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-elecs_fpath",
        required=True,
        help="The input BIDS `*electrodes.tsv` file with all "
             "the electrode contacts labeled and localized. "
             "A corresponding `*coordsystem.json` file must "
             "also accompany it.",
    )
    parser.add_argument(
        "-coord_frame",
        help="The output BIDS coordinate frame.",
        required=True,
    )
    parser.add_argument(
        "-output_elecs_fpath",
        help="The output BIDS datafile for electrodes in tsv format.",
        required=True,
    )
    parser.add_argument(
        "-bids_root",
        help="The root of the BIDS dataset.",
        required=True,
    )
    parser.add_argument(
        "-subjects_dir",
        help="The root of the BIDS dataset.",
        required=True,
    )

    args = parser.parse_args()

    # Extract arguments from parser
    elecs_fpath = args.elecs_fpath
    coord_frame = args.coord_frame
    output_elecs_fpath = args.output_elecs_fpath
    bids_root = args.bids_root
    subjects_dir = args.subjects_dir

    # run main script
    main(bids_root, elecs_fpath, output_elecs_fpath,
         to_coord_frame=coord_frame, subjects_dir=subjects_dir)

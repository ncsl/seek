import argparse
import sys
import json

sys.path.append("../../")

from neuroimg.base.objects.neuroimaging.fs_stats_objects import (
    FSCorticalStats,
    FSSegmentationStats,
)
from neuroimg.base.objects.neuroimaging.baseneuroimage import (
    RegionIndexMapping,
    ColorLut,
)


def create_mrtrixlut_volumes(fslut, mrtrixlut, regionmapping, fsvolumedict):
    volumes = {}

    for name in fsvolumedict.keys():
        # find name in FS LUT File
        src_lut_ind = [
            fslut.inds[ind]
            for ind, structname in enumerate(fslut.names)
            if name in structname.lower()
        ]

        if len(src_lut_ind) == 0:
            src_lut_ind = None
        else:
            src_lut_ind = src_lut_ind[0]

        # get corresponding index in Target LUT File
        trg_lut_ind = regionmapping.source_to_target(src_lut_ind)

        # create mapping from lut target index to the volume
        volumes[str(trg_lut_ind)] = int(fsvolumedict[name])

    return volumes


def create_volume_dict(fsobjs, ordering=[]):
    volumes = {}

    for j, obj in enumerate(fsobjs):
        if ordering:
            currorder = ordering[j]
        else:
            currorder = ""

        for idx, name in enumerate(obj.names):
            volumes[currorder + name] = obj.volumes[idx]

    return volumes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "lh_atlas_stats", help="The stats file produced for this atlas."
    )
    parser.add_argument(
        "rh_atlas_stats", help="The stats file produced for this atlas."
    )
    parser.add_argument(
        "subcort_atlas_stats",
        help="The stats file produced for the subcortical segmentation",
    )
    parser.add_argument(
        "wm_stats", help="The stats file produced for the white matter detected."
    )

    parser.add_argument("fs_lut", help="Brain parcellation lut for FS.")
    parser.add_argument(
        "fsl_lut", help="Brain parcellation look up table to map the correct labels."
    )
    parser.add_argument("outputfile", help="Output file path to save results.")
    args = parser.parse_args()

    # extract arguments from parser
    lh_atlas_stats_file = args.lh_atlas_stats
    rh_atlas_stats_file = args.rh_atlas_stats
    subcort_stats_file = args.subcort_atlas_stats
    wm_stats_file = args.wm_stats

    fslutfile = args.fs_lut
    mrtrixlutfile = args.fsl_lut
    outputvolumefile = args.outputfile

    # create table objects for lh/rh/subcort and wm stats files
    lhstats_df = FSCorticalStats(lh_atlas_stats_file)
    rhstats_df = FSCorticalStats(rh_atlas_stats_file)
    subcortstats_df = FSSegmentationStats(subcort_stats_file)
    wmstats_df = FSSegmentationStats(wm_stats_file)

    # create lh/rh/subcort volume dictionary
    fsvolumedict = create_volume_dict([lhstats_df, rhstats_df, subcortstats_df])

    # initialize LUT file objects
    fslut = ColorLut(fslutfile)
    mrtlut = ColorLut(mrtrixlutfile)
    regionmapping = RegionIndexMapping(fslutfile, mrtrixlutfile)

    # create volume dictionary for mrtrix labels
    mrtrixvolumedict = create_mrtrixlut_volumes(
        fslut, mrtlut, regionmapping, fsvolumedict
    )

    print(mrtrixvolumedict)
    # save the mrtrixvolumedictionary in mm^3
    with open(outputvolumefile, "w", encoding="utf8") as fp:
        json.dump(
            mrtrixvolumedict,
            fp,
            indent=4,
            sort_keys=True,
            separators=(",", ": "),
            ensure_ascii=False,
        )

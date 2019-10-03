import argparse
import subprocess
import scipy.io
from pprint import pprint


def read_label_coords(elecfile):
    labels = []
    labelsxyz = []

    print("Reading ", elecfile)

    with open(elecfile, "r") as f:
        for i, _row in enumerate(f):
            if i == 0:
                hdr = _row.split(",")
                continue

            row = _row.split(",")
            labels.append(row[0])
            labelsxyz.append([float(x) for x in row[1:]])

    return labels, labelsxyz


def get_best_nghbr_match(all_tal_label):
    maxweight = 0

    for i in range(0, len(all_tal_label)):
        label = all_tal_label[i]

        weight, labeling = label.split(":")
        hemi, lobe, tal_label, matter, brodmanarea = labeling.split(",")

        if tal_label == "*":
            continue

        print(f"Current {weight} and {maxweight}")
        if float(weight) > maxweight:
            bestfound = (hemi, lobe, tal_label, matter, brodmanarea)
            # print(bestfound)
            maxweight = float(weight)

    if maxweight == 0:
        return maxweight, "*", "*", "*", "*", "*"

    # extract the meaning from best found labeling
    hemi, lobe, tal_label, matter, brodmanarea = bestfound
    return maxweight, hemi, lobe, tal_label, matter, brodmanarea


def nearest_gm_search(labels, labelsxyz):
    """
    Performs nearest Gray-Matter search for atlas label in Talairach space using the xyz coordinates (in mm)
    of contacts passed in.

    Essentially subprocess runs: "java -cp ./talairach.jar org.talairach.PointToTD 1, <points>;"

    Ref: http://www.talairach.org/manual.html

    :param labels:
    :param labelsxyz:
    :return:
    """
    assert len(labels) == len(labelsxyz)

    # initialize list to store labels
    atlaslabels = []
    hemilabels = []
    lobelabels = []
    matterlabels = []
    brodmanlabels = []

    # go through each coord one by one and run subprocess to extract talairach label
    for i in range(len(labels)):
        # get the x,y,z coordinates
        labelx, labely, labelz = labelsxyz[i]
        print(labelsxyz[i])
        # run subprocess to extract closest gm point in talairach atlas
        cp = subprocess.run(
            "java -cp ./talairach.jar org.talairach.PointToTD 3:9, %s, %s, %s"
            % (labelx, labely, labelz),
            shell=True,
            stdout=subprocess.PIPE,
        )

        # label
        all_tal_labels = cp.stdout.decode("ascii").strip().split("\n")

        all_tal_labels = [x.strip() for x in all_tal_labels][2:]
        # print(len(all_tal_labels))
        pprint(all_tal_labels)
        # pprint(f"Tal label: \n {all_tal_labels} \nfor {labels[i]}")

        weight, hemi, lobe, tal_label, matter, brodmanarea = get_best_nghbr_match(
            all_tal_labels
        )

        print(f"Found: {weight}, {hemi}, {lobe}, {tal_label}, {matter}, {brodmanarea}")

        hemilabels.append(hemi)
        matterlabels.append(matter)
        lobelabels.append(lobe)
        atlaslabels.append(tal_label)
        brodmanlabels.append(brodmanarea)

    return atlaslabels


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("txt_talairach_file")
    parser.add_argument("labeled_txt_talairach_file")
    args = parser.parse_args()

    # extract arguments from parser
    txt_talairach_file = args.txt_talairach_file
    labeled_txt_talairach_file = args.labeled_txt_talairach_file

    # read in electrodes file
    labels, labelsxyz = read_label_coords(txt_talairach_file)

    print("Found labels: ", labels, "\n\n")
    # print(labelsxyz)
    # run coordinate transformation labeling using talairach client
    atlaslabels = nearest_gm_search(labels, labelsxyz)

    mdict = {"eleclabels": labels, "elecmatrix": labelsxyz, "anatomy": atlaslabels}

    scipy.io.savemat(labeled_txt_talairach_file, mdict=mdict)

# -*- coding: utf-8 -*-

# Data structure manipulations and conversions

import bz2
import contextlib
import json
import pickle
import re
from collections import OrderedDict
from copy import deepcopy
from datetime import date, datetime

import numpy as np
import scipy.io


class MatReader:
    """
    Object to read mat files into a nested dictionary if need be.
    Helps keep strucutre from matlab similar to what is used in python.
    """

    def __init__(self, filename=None):
        self.filename = filename

    def loadmat(self, filename):
        """
        this function should be called instead of direct spio.loadmat
        as it cures the problem of not properly recovering python dictionaries
        from mat files. It calls the function check keys to cure all entries
        which are still mat-objects
        """
        data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
        return self._check_keys(data)

    def _check_keys(self, dict):
        """
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in dict:
            if isinstance(dict[key], scipy.io.matlab.mio5_params.mat_struct):
                dict[key] = self._todict(dict[key])
        return dict

    def _todict(self, matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        dict = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
                dict[strg] = self._todict(elem)
            elif isinstance(elem, np.ndarray):
                dict[strg] = self._tolist(elem)
            else:
                dict[strg] = elem
        return dict

    def _tolist(self, ndarray):
        """
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, scipy.io.matlab.mio5_params.mat_struct):
                elem_list.append(self._todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(self._tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list

    def convertMatToJSON(self, matData, fileName):
        jsonData = {}

        for key in matData.keys():
            if (type(matData[key])) is np.ndarray:
                serializedData = pickle.dumps(
                    matData[key], protocol=0
                )  # protocol 0 is printable ASCII
                jsonData[key] = serializedData
            else:
                jsonData[key] = matData[key]

        with contextlib.closing(bz2.BZ2File(fileName, "wb")) as f:
            json.dump(jsonData, f)


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(
                obj,
                (
                        np.int_,
                        np.intc,
                        np.intp,
                        np.int8,
                        np.int16,
                        np.int32,
                        np.int64,
                        np.uint8,
                        np.uint16,
                        np.uint32,
                        np.uint64,
                ),
        ):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):  # This is the fix
            return obj.tolist()
        elif isinstance(obj, (datetime, date)):
            return obj.isoformat()
        return json.JSONEncoder.default(self, obj)


def vector2scalar(x):
    if not (isinstance(x, np.ndarray)):
        return x
    else:
        y = np.squeeze(x)
    if all(y.squeeze() == y[0]):
        return y[0]
    else:
        return reg_dict(x)


def list_of_strings_to_string(lstr, sep=","):
    result_str = lstr[0]
    for s in lstr[1:]:
        result_str += sep + s
    return result_str


def dict_str(d):
    s = "{"
    for key, value in d.items():
        s += "\n" + key + ": " + str(value)
    s += "}"
    return s


def split_string_text_numbers(ls):
    items = []
    for s in ensure_list(ls):
        match = re.findall("(\d+|\D+)", s)
        if match:
            items.append(tuple(match[:2]))
    return items


def construct_import_path(path, package="tvb_epilepsy"):
    path = path.split(".py")[0]
    start = path.find(package)
    return path[start:].replace("/", ".")


def formal_repr(instance, attr_dict, sort_dict_flag=False):
    """ A formal string representation for an object.
    :param attr_dict: dictionary attribute_name: attribute_value
    :param instance:  Instance to read class name from it
    """
    class_name = instance.__class__.__name__
    formal = class_name + "{"
    if sort_dict_flag:
        attr_dict = sort_dict(attr_dict)
    for key, val in attr_dict.items():
        if isinstance(val, dict):
            formal += "\n" + key + "=["
            for key2, val2 in val.items():
                formal += "\n" + str(key2) + " = " + str(val2)
            formal += "]"
        else:
            formal += "\n" + str(key) + " = " + str(val)
    return formal + "}"


def obj_to_dict(obj):
    """
    :param obj: Python object to introspect
    :return: dictionary after recursively taking obj fields and their values
    """
    if obj is None:
        return obj
    if isinstance(obj, (str, int, float)):
        return obj
    if isinstance(obj, (np.float32,)):
        return float(obj)
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    if isinstance(obj, list):
        ret = []
        for val in obj:
            ret.append(obj_to_dict(val))
        return ret
    ret = {}
    for key in obj.__dict__:
        val = getattr(obj, key, None)
        ret[key] = obj_to_dict(val)
    return ret


def reg_dict(x, lbl=None, sort=None):
    """
    :x: a list or np vector
    :lbl: a list or np vector of labels
    :return: dictionary
    """
    if not (isinstance(x, (str, int, float, list, np.ndarray))):
        return x
    else:
        if not (isinstance(x, list)):
            x = np.squeeze(x)
        x_no = len(x)
        if not (isinstance(lbl, (list, np.ndarray))):
            lbl = np.repeat("", x_no)
        else:
            lbl = np.squeeze(lbl)
        labels_no = len(lbl)
        total_no = min(labels_no, x_no)
        if x_no <= labels_no:
            if sort == "ascend":
                ind = np.argsort(x).tolist()
            elif sort == "descend":
                ind = np.argsort(x)
                ind = ind[::-1].tolist()
            else:
                ind = range(x_no)
        else:
            ind = range(total_no)
        d = OrderedDict()
        for i in ind:
            d[str(i) + "." + str(lbl[i])] = x[i]
        if labels_no > total_no:
            ind_lbl = np.delete(np.array(range(labels_no)), ind).tolist()
            for i in ind_lbl:
                d[str(i) + "." + str(lbl[i])] = None
        if x_no > total_no:
            ind_x = np.delete(np.array(range(x_no)), ind).tolist()
            for i in ind_x:
                d[str(i) + "."] = x[i]
        return d


def sort_dict(d):
    return OrderedDict(sorted(d.items(), key=lambda t: t[0]))


def dicts_of_lists(dictionary, n=1):
    for key, value in dictionary.items():
        dictionary[key] = ensure_list(dictionary[key])
        if len(dictionary[key]) == 1 and n > 1:
            dictionary[key] = dictionary[key] * n
    return dictionary


def iterable_to_dict(obj):
    d = OrderedDict()
    for ind, value in enumerate(obj):
        d["%02d" % ind] = value
    return d


def dict_to_list_or_tuple(dictionary, output_obj="list"):
    dictionary = sort_dict(dictionary)
    output = dictionary.values()
    if output_obj == "tuple":
        output = tuple(output)
    return output


def list_of_dicts_to_dicts_of_ndarrays(lst, shape=None):
    d = dict(zip(lst[0], zip(*list([d.values() for d in lst]))))
    if isinstance(shape, tuple):
        for key, val in d.items():
            d[key] = np.reshape(np.stack(d[key]), shape)
    else:
        for key, val in d.items():
            d[key] = np.squeeze(np.stack(d[key]))
    return d


def arrays_of_dicts_to_dicts_of_ndarrays(arr):
    lst = arr.flatten().tolist()
    d = list_of_dicts_to_dicts_of_ndarrays(lst)
    for key, val in d.items():
        d[key] = np.reshape(d[key], arr.shape)
    return d


def dicts_of_lists_to_lists_of_dicts(dictionary):
    return [dict(zip(dictionary, t)) for t in zip(*dictionary.values())]


def ensure_string(arg):
    if not (isinstance(arg, str)):
        if arg is None:
            return ""
        else:
            return ensure_list(arg)[0]
    else:
        return arg


def ensure_list(arg):
    if not (isinstance(arg, list)):
        try:  # if iterable
            if isinstance(arg, (str, dict)):
                arg = [arg]
            else:
                arg = list(arg)
        except BaseException:  # if not iterable
            arg = [arg]
    return arg


def ensure_string(arg):
    if not (isinstance(arg, str)):
        if arg is None:
            return ""
        else:
            return ensure_list(arg)[0]
    else:
        return arg


def set_list_item_by_reference_safely(ind, item, lst):
    while ind >= len(lst):
        lst.append(None)
    lst.__setitem__(ind, item)


def get_list_or_tuple_item_safely(obj, key):
    try:
        return obj[int(key)]
    except BaseException:
        return None


def linear_index_to_coordinate_tuples(linear_index, shape):
    if len(linear_index) > 0:
        coordinates_tuple = np.unravel_index(linear_index, shape)
        return zip(*[ca.flatten().tolist() for ca in coordinates_tuple])
    else:
        return []


def labels_to_inds(labels, lbls):
    idx = []
    lbls = ensure_list(lbls)
    for i, label in enumerate(labels):
        if label in lbls:
            idx.append(i)
    return np.unique(idx)


def generate_region_labels(n_regions, labels=[], str=". ", numbering=True):
    if len(labels) == n_regions:
        if numbering:
            return np.array(
                [
                    str.join(["%d", "%s"]) % tuple(l)
                    for l in zip(range(n_regions), labels)
                ]
            )
        else:
            return labels
    else:
        return np.array(["%d" % l for l in range(n_regions)])


def monopolar_to_bipolar(labels, indices=None, data=None):
    if indices is None:
        indices = range(len(labels))
    bipolar_lbls = []
    bipolar_inds = [[], []]
    for ind in range(len(indices) - 1):
        iS1 = indices[ind]
        iS2 = indices[ind + 1]
        if (labels[iS1][0] == labels[iS2][0]) and int(
                re.findall(r"\d+", labels[iS1])[0]
        ) == int(re.findall(r"\d+", labels[iS2])[0]) - 1:
            bipolar_lbls.append(labels[iS1] + "-" + labels[iS2])
            bipolar_inds[0].append(iS1)
            bipolar_inds[1].append(iS2)
    if isinstance(data, np.ndarray):
        data = data[bipolar_inds[0]] - data[bipolar_inds[1]]
        return bipolar_lbls, bipolar_inds, data
    else:
        return bipolar_lbls, bipolar_inds


def shape_to_size(shape):
    shape = np.array(shape)
    shape = shape[shape > 0]
    return np.int(np.max([shape.prod(), 1]))


def shape_to_ndim(shape, squeeze=False):
    if squeeze:
        shape = filter(lambda x: not (np.any(np.in1d(x, [0, 1]))), list(shape))
    return len(shape)


def linspace_broadcast(start, stop, num_steps, maxdims=3):
    x_star = np.linspace(0, 1, num_steps)
    dims = 0
    x = None
    while x is None and dims < maxdims:
        try:
            x = x_star[:, None] * (stop - start) + start
        except BaseException:
            x_star = x_star[:, np.newaxis]
            dims = dims + 1
    return x


def squeeze_array_to_scalar(arr):
    arr = np.array(arr)
    if arr.size == 1:
        return arr
    elif np.all(arr == arr[0]):
        return arr[0]
    else:
        return arr


def make_float(x, precision="64"):
    if isinstance(x, np.ndarray):
        if isequal_string(precision, "64"):
            return x.astype(np.float64)
        elif isequal_string(precision, "32"):
            return x.astype(np.float32)
        else:
            return x.astype(np.float)
    else:
        if isequal_string(precision, "64"):
            return np.float64(x)
        elif isequal_string(precision, "32"):
            np.float32(x)
        else:
            return np.float(x)


def make_int(x, precision="64"):
    if isinstance(x, np.ndarray):
        if isequal_string(precision, "64"):
            return x.astype(np.int64)
        elif isequal_string(precision, "32"):
            return x.astype(np.int32)
        else:
            return x.astype(np.int)
    else:
        if isequal_string(precision, "64"):
            return np.int64(x)
        elif isequal_string(precision, "32"):
            np.int32(x)
        else:
            return np.int(x)


def copy_object_attributes(
        obj1, obj2, attr1, attr2=None, deep_copy=False, check_none=False
):
    attr1 = ensure_list(attr1)
    if attr2 is None:
        attr2 = attr1
    else:
        attr2 = ensure_list(attr2)
    if deep_copy:

        def fcopy(a1, a2):
            return setattr(obj2, a2, deepcopy(getattr(obj1, a1)))

    else:

        def fcopy(a1, a2):
            return setattr(obj2, a2, getattr(obj1, a1))

    if check_none:
        for a1, a2 in zip(attr1, attr2):
            if getattr(obj2, a2) is None:
                fcopy(a1, a2)
    else:
        for a1, a2 in zip(attr1, attr2):
            fcopy(a1, a2)
    return obj2

import datetime
import json
from datetime import date

import numpy as np
from enum import Enum

try:
    to_unicode = unicode
except NameError:
    to_unicode = str


class ResultTypes(Enum):
    mvar = 'mvar'
    pert = 'pert'
    stft = 'stft'


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):  # This is the fix
            return obj.tolist()
        elif isinstance(obj, (datetime, date)):
            return obj.isoformat()
        return json.JSONEncoder.default(self, obj)


# possible reference schemes for scalp recording data - note that there is no bipolar
SCALP_REFERENCE_CONFIGS = [
    'monopolar',
    'common_avg'
]

IEEG_REFERENCE_CONFIGS = [
    'monopolar',
    'bipolar',
    'common_avg'
]

# POSSIBLE DATA TYPES IN RECORDING DATA
ICTAL_TYPES = ['ictal', 'sz', 'seiz', 'seizure']
INTERICTAL_TYPES = ['interictal', 'ii', 'aslp', 'aw']

RAWMODALITIES = ['ieeg', 'seeg', 'scalp', 'ecog']

# POSSIBLE DATA TYPES IN RESULTS MODELS
MODELDATATYPES = ['ltv', 'pert', 'frag', 'freq']

# frequency band
DALPHA = [0, 15]
BETA = [15, 30]
GAMMA = [30, 90]
HIGH = [90, 200]

class Freqbands(Enum):
    DALPHA = DALPHA
    BETA = BETA
    GAMMA = GAMMA
    HIGH = HIGH
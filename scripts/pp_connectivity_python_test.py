import scipy.io
import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.datasets import sample

from mne_connectivity import spectral_connectivity_epochs



mat = scipy.io.loadmat('file.mat')



con_wpli = spectral_connectivity_epochs(
    epochs,
    method="wpli",
    mode="multitaper",
    sfreq=sfreq,
    fmin=fmin,
    fmax=fmax,
    faverage=True,
    tmin=tmin,
    mt_adaptive=False,
    n_jobs=1,
)
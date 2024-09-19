import scipy.io
import matplotlib.pyplot as plt
import mne
import numpy as np
from mne_connectivity import spectral_connectivity_epochs





epochs = scipy.io.loadmat('C:/Users/timdr/OneDrive/Uni_Oldenburg/3_Semester/Module/Pratical_Project/Analysis/data/sandbox/testdata.mat')

print(type(epochs))




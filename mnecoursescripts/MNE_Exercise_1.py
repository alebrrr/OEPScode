# Now that you have all the nesseccary tools, you can apply them to actual
# data. The aim of this exercise is that you reproduce what we did together
# using a new dataset - i will tell you the steps you have to do but you
# need to implement them by yourself

# STEP1: Load the data
# this should be easy you already did that. You can choose any subject.

# STEP2: filter the data, you have to decide on the cutoff frequencies

# STEP3: ICA, do an independent component analysis and remove blinks.
# (If you don't manage to do that you can skip this step)

# STEP4: Epoch the data. Now you have to decide what your epoch parameters are
# going to be. Paramters to consider are:
# -start and stop (tmin/tmax)
# -baseline
# -rejection threshold
# Apparently for some subjects the some event appears twice in the list
# which causes an error in the epochs function. You can set the parameter
# event_repeated to "drop" which means that the duplet is being dropped.

# STEP5: Analyse ERPs. Plot the evoked response for the auditory stimulus,
# visual stimulus, auditory omission and visual omission. Also plot the
# topography. Compare them. How could you quantify the evoked responses?

import numpy as np
import mne
list_of_subjects = ["001", "0002", "003"]
#reference_ica = mne.preprocessing.read_ica("")

paths = []
for subject in list_of_subjects:
    path = "/home/bialas/projects/eeg_course/ds002218-download/sub-%s/eeg/sub-%s_task-Experiment_eeg.set" % (
        subject, subject)
    paths.append(path)

    # raw =
    # raw.filter
    # ica.
    # corrmap
    # ica.apply()
    # autoreject
    # epochs.save("sub-%s/eeg/sub-%s_task-Experiment_eeg.set" % (
    #    subject, subject))

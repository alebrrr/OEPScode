**2021 ale b code for OEPS data to MNE pipeline**

preprocessing-pipeline contains all functions necessary to go from OEPS binaries to MNE objects.

When recording OEPS creates a folder with input subject name and timestamp on OEPS computer. 
This folder contains all recordings for one continuous data streaming session (one recording= one labplatform experiment, the experiment starts and stops recording with sustained trigger). 

in OEPS'structure "experiment" is one data streaming session, with multiple recordings, which is confusing. in the code I refer to OEPS "experiments" as "sessions" but I do not change folder naming.

prior to any processing the subject named OEPS folder with all its subfolders is moved to the storage (I am putting them in alessandro/mea). h5 file on tdt computer is moved in corresponding "recording" folder.
prova is a binary of downsampled data. OpenEphys.py one of the original OEPS scrypt for reading data. run_pipeline is scraps, code I used to make plots and calculate snr in a pinch.

channelsH32M contains relative channel positions (sphere with r 1). channelsH32Mabs contains positions in millimiters.

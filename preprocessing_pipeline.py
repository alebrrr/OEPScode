from mne import set_eeg_reference, events_from_annotations
from mne.epochs import Epochs
# from autoreject import AutoReject, Ransac
import numpy as np
from mne.preprocessing.ica import ICA, corrmap, read_ica
import os
from pathlib import Path
from matplotlib import pyplot as plt, patches
import mne


def loadOep(folder="Z:\\Alessandro_Braga\\MEA data february\\feb_n1_2021-02-26_16-35-13_m0002\\Record Node 101",session=1, newfreq=1000, recording=1, mainevent=8):
    """
    all these functions rely on filenames and folders retaining the structure set by oeps. h5 files of a certain recording-experiment must be put in the corresponding recordingN (e.g. recording5) folder.

    This function lods OEPS data. folder structure needs to be the original OEPS one. Mainevent is key for trial presentation indipendent from  experiment type
    

    Parameters
    ----------
    folder : str, folder of recording session (sereis of recordings) up to subfolder Record Node 101".
    session :int, which session (if there are more than one  session for animal, that is, if oeps stream is turned off and back on)
    newfreq : int , freq to downsample data to
        
    recording : int: which of various recordings during one data streaming session. each recording is one experiment
        DESCRIPTION. The default is 10.
    mainevent : int
        Mainevent is key for trial presentation indipendent from  experiment type. default 8

    Returns (#TODO CHECK WHICH INTS ARE ACTUALLY FLOATS)
    -------
    indexo : array or []
        stimulus index if mmn from h5..
    filenamex : filename for recording being analised
    soa : float. TODO get to int
        soa from h5. if h5 damaged it is calculated from n trial events and recording lenght 
    stimdur float. TODO get to int
        stimulus duration  from h5. if h5 damaged it is calculated from n trial events and recording lenght 
    mean : array
        mean=averaged epoched traces.
    datar : array
        downsampled data with channels averaged together to reduce channel count to 6. does not work well, pooled channels are noisy.
    data
    : array
        downsampled data.
    triggers : array 
        trigger events index for timestamps in teimestams.
    timestamps : array
        timestamps , need zeroing.
    sfreq : int
               sfrequ is sampling freq should always be newfreq at this pOint
.
    oldfreq : int
        original freq.

    """
    import numpy as np
    from ast import literal_eval
    from glob import glob
    global _out_folder
    global front_left
    global media_left
    global poste_left
    global front_right
    global media_right
    global poste_right
    global cent
    global omiss
    global poste_cent
    omiss = False
    indexo = []  # it stays [] if all trials are identical and there is no index
    front_left = ['B6', 'B5', 'C6', 'C5']  # for future use
    channelmap = np.array(
        [[16, 18, 15, 17], [14, 12, 13, 11], [10, 0, 9, 7], [20, 29, 19, 22], [6, 3, 5, 2], [23, 27, 24, 26],
         [8, 1, 21, 28], [28, 1, 4, 25]])  # see below

    front_right = ['D6', 'D5', 'E6', 'E5']
    front_righti = [14, 12, 13, 11]
    media_right = ['D4', 'F3', 'E4', 'E3']
    media_righti = [10, 0, 9, 7]
    media_left = ['B4', 'A3', 'C4', 'B3']
    media_lefti = [20, 29, 19, 22]
    poste_right = ['F1', 'F2', 'E1', 'E2']
    poste_righti = [6, 3, 5, 2]
    poste_left = ['A1', 'B2', 'B1', 'A2']
    poste_lefti = [23, 27, 24, 26]
    cent = ['D3', 'C2', 'C3', 'D2']
    centi = [8, 1, 21, 28]
    poste_cent = ['D1', 'C2', 'C1', 'D2']
    poste_centi = [28, 1, 4, 25]
    _out_folder = 'C:\\Users\PC\Desktop\mea_anal'

    Folder = folder  # folder where the recordings are
    session = session  # usually 1. it is the overall recording session
    Recording = recording  # various recordings per session
    Processor = None
    Files = sorted(glob(Folder + '/**/*.dat', recursive=True))  # get all data files in folder
    InfoFiles = sorted(glob(Folder + '/*/*/structure.oebin'))

    Data, Rate = {}, {}
    for F, File in enumerate(Files):  # go through, select the recording  from input.

        Exp, Rec, _, Proc = File.split('\\')[-5:-1]
        Exp = str(int(Exp[10:]) - 1)
        Rec = str(int(Rec[9:]) - 1)
        Proc = Proc.split('.')[0].split('-')[-1]
        if '_' in Proc: Proc = Proc.split('_')[0]

        if Proc not in Data.keys(): Data[Proc], Rate[Proc] = {}, {}

        if session:
            if int(Exp) != session - 1: continue

        if Recording:
            if int(Rec) != Recording - 1: continue

        if Processor:
            if Proc != Processor: continue

        print('Loading recording', int(Rec) + 1, '...')
        if Exp not in Data[Proc]: Data[Proc][Exp] = {}
        Data[Proc][Exp][Rec] = np.memmap(File, dtype='int16', mode='c')

        Info = literal_eval(open(InfoFiles[F]).read())
        ProcIndex = [Info['continuous'].index(_) for _ in Info['continuous']
                     if str(_['recorded_processor_id']) == '101'][0]

        ChNo = Info['continuous'][ProcIndex]['num_channels']
        if Data[Proc][Exp][Rec].shape[0] % ChNo:
            print('Rec', Rec, 'is broken')
            del (Data[Proc][Exp][Rec])
            continue

        SamplesPerCh = Data[Proc][Exp][Rec].shape[0] // ChNo
        Data[Proc][Exp][Rec] = Data[Proc][Exp][Rec].reshape((SamplesPerCh, ChNo))
        Rate[Proc][Exp] = Info['continuous'][ProcIndex]['sample_rate']
        f = Rate[Proc][Exp]
        oldfreq = f
        FolderEv = Folder + '\\experiment' + str(session) + '\\recording' + str(
            recording) + '\\events\Rhythm_FPGA-100.0\\TTL_1'  # get folder with events
        FilesEv = sorted(glob(FolderEv + '\\*.npy', recursive=True))  # get events
        ttll = np.load(FilesEv[0])

        ttl = np.argwhere(ttll == mainevent)  # get trial event.
        tstp = np.load(FilesEv[3])  # load timestamps
        data = Data[Proc][Exp][str(recording - 1)]
        chunk = int((len(Data[Proc][Exp][str(recording - 1)]) / len(ttl)) / (int(f / newfreq)))  # size of trial in datapoints

        filenamex = File  # to set aside the filename of the actual recording we are working with after we exit the loop through all FIles for that session

        u = []

    Data = []

    myfile = filenamex[:-3] + str(
        newfreq) + '_dwnsmpl.npy'  # this is to average together channels as indicated in the top list (channelmap). looked weird when i actually used these data TODO look into averaging together neighboring channels to go from 30 to 8

    from pathlib import Path

    my_file = Path(myfile)
    if my_file.is_file():
        print('load downsampled file')

        data = np.load(myfile)
        f = newfreq
    else:
        print('saving downsampled file')  # if frequency is not newfreq hz data are downsampled to newfreq
        if f != newfreq:
            data = data[::int(f / newfreq), :32]  # downsample from 30k
            f = newfreq
        data = data[:, :32]  # cut dac channels

        data = np.delete(data, [7, 24], 1)  # get rid of non-channels

        data = (
                    data * 0.195 / 1000000)  # convert bits to microvolt and microvolt to volt. "The raw data are saved as signed 16-bit integers, in the range -32768 to 32767. They don’t have a unit. To convert to microvolts, just  multiply by 0.195, as you noted. This scales the data into the range ±6.390 mV, with 0.195 µV resolution (Intan chips have a ±5 mV input range).""
        # THIS GETS DATA IN ACTUAL MEMORY; TAKES TIME
        data = data - data.mean(0)  # baseline
        np.save(myfile, data)

    h5file = sorted(glob(Folder + '/experiment' + str(session) + '/recording' + str(recording) + '/**/*.h5',
                         recursive=True))  # this is the file from the tdt side. they are put by hand in the recordingn folder after data collection
    import tables
    try:
        ff = tables.open_file(h5file[0])
        soa = ((
                    ff.root._v_attrs.Exp_Settings_SOA / ff.root._v_attrs.Exp_Settings_sampling_freq))  # SOA for filenames. rounding issue here
        stimdur = ((
                    ff.root._v_attrs.Exp_Settings_stim_dur / ff.root._v_attrs.Exp_Settings_sampling_freq))  # stim duration for filenames
        if h5file[0][-14] == 'N':  # checks if the recording is a MMN one, in which case gets the index.
            omiss = True
            indexo = Omitter(h5file[
                                 0])  # get index of omissions from recorded audio signal. this allows to check for discrepancies from index (there should be none if the computer is not used during recording)
    except Exception:  # in case the h5 file is damaged or does not exist.
        soa = ((len(data) / (newfreq ** 2)))
        stimdur = 0.109
    tstp = tstp - tstp[0]
    index = tstp[ttl] / (oldfreq / newfreq)  # downsample timestamps
    if index[-1] + chunk > len(data):  # hack  if recording happens to be a datapoint too short last trial is deleted.
        index = np.delete(index, len(index) - 1)
        deletable = np.argwhere((tstp - tstp[0]) / int(oldfreq / newfreq) + chunk > len(data))
        tstp = np.delete(tstp, deletable)
        ttll = np.delete(ttll, deletable)
    v = []  # these following loops are to average epochs , the purpose being having quick way to look at erp after recording session

    for iz in range(0, len(data[1, :])):
        n = []

        for i in range(0, len(index)):
            n.append(data[int(index[i]):int(index[i]) + chunk, iz])

        v.append(sum(n) / len(n) - sum(sum(n)) / (len(n) * len(sum(n))))

    mean = v  # this is the epoched and averaged data
    triggers = ttll
    sfreq = f
    timestamps = tstp
    return indexo, filenamex, soa, stimdur, mean, data, triggers, timestamps, sfreq, oldfreq


# indexo: stimulus index if mmn from h5. filename. stimulus onset async from h5. stimulus duration from h5. mean=averaged epoched traces, data = downsampled data, datar= downsampled data iwth channels averaged in groups (left front ecc). triggers= trigger events index for timestamps in teimestams, sfrequ is sampling freq should always be newfreq at this pint


def Omitter(h5file):
    """
    this method goes through the audio recording and checks which trials are silent. this is a temporary solution for omission MMN as in two recordings the index was not followed (intensity switching logic is time dependent and can fail. Should be ok now)
    """
    omiss: True

    filename = h5file

    import tables

    import numpy as np
    import numpy

    ff = tables.open_file(filename)
    leftn1 = []
    for group in ff.walk_groups("/"):
        for array in ff.list_nodes(group, classname='Array'):
            leftn1.append(array[:, 1])
    leftn1 = np.array(leftn1) * 50000
    leftarr = np.zeros((ff.root._v_attrs.Exp_Settings_reps, int(ff.root._v_attrs.Exp_Settings_SOA / 10)))
    for i in range(0, ff.root._v_attrs.Exp_Settings_reps):
        leftarr[i] = numpy.array(leftn1)[0, (i * int(ff.root._v_attrs.Exp_Settings_SOA / 10)): (i + 1) * int(
            ff.root._v_attrs.Exp_Settings_SOA / 10)]
    dd = np.max(leftarr[:, 375:625], 1)
    uu = dd > 100
    return uu


def scaler(my_diction, factor=50):
    total = 0

    for j in my_diction:
        my_diction[j] = (my_diction[j]) / factor
    return my_diction


def MNEify(indexo, datar, timestamps, triggers, sfreq, oldfreq, trialevent=8, scale_factorX=50, scale_factorDIV=50):
    """
            this method gets the data into an Mne raw object


    Parameters
    ----------
    indexo :array, index of different events. . is [] if all trials are the same
    datar : array, data
        .
    timestamps :array, timestamps obtained from file as generated by Oeps
    triggers : corresponding events from timestamps
    sfreq : data sampling frequency, int
    oldfreq : original data sampling frequency. 
    trialevent : what int corresponds to trial event in triggers
    scale_factorX:factor to multiply nas and ear points for digmontage. (now expressed in meters, with actual values in the order of millimiters)
    scale_factorDIV:factor to input to scaler for the channel positions. (now expressed in mm) 
    Returns
    -------
    raw objects and event array to make epochs
    """
    # TODO. RIGHT NOW i AM DIVIDING THE ch position values IN millimiters BY 50, and also  multiplying the nas +ear values in meters by 50, in order to get a sensible topomap.   counterintuitive.
    import mne
    datad = np.transpose(datar[:, :len(datar[1, :])])  # transpose raw data as per mne requirment
    timstp = timestamps - timestamps[0]  # zero timestamps to start of recording
    eventstsp = timstp[np.argwhere(triggers == trialevent)]  # select timestamps of each trial event
    ggg = np.zeros((len(eventstsp), 1))
    dd = np.ones((len(eventstsp), 1))
    if np.array(indexo).any():  # fix discrepancies between h5file index and timestamps
        dd = indexo.astype(int)
        dd = np.reshape(dd, (len(dd), 1))
        if len(dd) > len(eventstsp):
            dd = dd[:len(eventstsp)]
        elif len(dd) < len(eventstsp):
            eventstsp = eventstsp[:len(dd)]
            ggg = ggg[:len(dd)]

    #     tim=tim-tim[1843]#3000 datapoint at 48828.125hz recorded at 30000hz

    # unfilt=np.delete(unfilt,[7,24],0)
    # filtt=np.delete(filtt,[7,24],0)
    if len(datad) == 30:  # channel names for all 30
        import pandas as pd

        df = pd.read_csv('channelsH32Mabs.txt')
        ch_names = df.name.to_list()

        pos = df[['z', 'x', 'y']].values
        dig_ch_pos = dict(zip(ch_names, pos))
        scaler(dig_ch_pos, scale_factorDIV)
        nasion = [-0.007 * scale_factorX, 0 * scale_factorX, -0.001 * scale_factorX]
        lpa = [0.000 * scale_factorX, -0.005 * scale_factorX, -0.001 * scale_factorX]
        rpa = [0.000 * scale_factorX, 0.005 * scale_factorX, -0.001 * scale_factorX]
        montage = mne.channels.make_dig_montage(ch_pos=dig_ch_pos, nasion=nasion, lpa=lpa, rpa=rpa)

        ch_types = ['eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg',
                    'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg',
                    'eeg', 'eeg']
    else:  # channel names for pooled channels
        ch_names = ['front_left', 'front_right', 'cent_right', 'cent_left', 'post_right', 'post_left', 'central',
                    'post_cent']
        ch_types = ['eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg']

    events = np.concatenate((eventstsp / (oldfreq / sfreq) + int(np.ceil((310 / 4882.8125) * sfreq)), ggg, dd),
                            1)  # here we get events in shape required from mne (1st column timestamps, second zeros, third ones)
    events = events.astype(int)
    info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
    raw = mne.io.RawArray(datad, info)
    raw = raw.set_montage(montage)
    events = events.astype(int)
    return raw, events


def find_peaks(data, min_dist=1, thresh=0.3, degree=None):
    """
    Find peaks in the data, optionally apply a baseline before.
    min_dist: minimum distance between peaks, the biggest peak is preferred.
    thresh: float between 0.0 and 1.0, normalized detection threshold
    degree: degree of the polynomial that will estimate the data baseline.
    if None no baseline detection is done
    """
    import peakutils
    if degree is not None:
        base = peakutils.baseline(data, degree)
    else:
        base = np.zeros(len(data))
    peak_idx = peakutils.indexes(data - base, thres=thresh, min_dist=min_dist)
    if peak_idx.size == 0:
        peak_idx = np.argwhere(data == max(data))
    return peak_idx, base


def peak_minmax(trace, pz=[60, 140], nu=[160, 240]):
    """
        see findpeax above, this one to get adjacent positive and negative peaks to measure amplitude. did not test extensively yet


    Parameters
    ----------
    trace : array,evoked response
    pz :interval for positive peak, list of int
    nu : interval for negative peak, list of int

    Returns
    -------
    pospeak : 
        .
    negpeak :int, neg peak voltage
    pospeakindex : 
    negpeakindex : ints, peaks position from first value of scanned interval

    """

    bb, roe = find_peaks(trace[pz[0]:pz[1]], min_dist=1, thresh=0.3, degree=None)
    tt, roe = find_peaks(0 - trace[nu[0]:nu[1]], min_dist=1, thresh=0.3, degree=None)
    pospeakindex = pz[0]+(bb[max(trace[pz[0]:pz[1]][bb])==trace[pz[0]:pz[1]][bb]])
    pospeak = trace[pospeakindex]
    negpeakindex = nu[0] + (tt[max(0-trace[nu[0]:nu[1]][tt])==0-trace[nu[0]:nu[1]][tt]])
    negpeak = trace[negpeakindex]
    ampl=max([pospeak,negpeak])-min([pospeak,negpeak])

    return ampl, pospeakindex


def noise_rms(epochedt, chan, tr):
    """
    Calculate the noise that remains after averaging in the ERP following
    the procedure from Schimmel(1967): invert every other trial then average.
    Afterwards calculate the root mean square across the whole EPR.
    
    parameters
    --------
        epochedt= copy of array obtained form mne epochs . epochs['stim_standard']) (if epochs from different events otherwise  np.array(epochs) + epoched.copy())
        chan= int, which channel is being worked with
        tr= int, up to wich trial is rms being calculated
    returns
    --------
    rms for a certain channel evoked response over a certain amount of trials
    
    """

    for i in range(tr):
        if not i % 2:
            epochedt[i, :, :] = -epochedt[i, :, :]
    evokedchannel = np.array(epochedt)[0:tr, chan, :].mean(0)
    rms = np.sqrt(np.mean(evokedchannel ** 2))

    return rms


def noise_inverter(epochedt):
    """
    see noise rms but without actual calculation
    
    Parameters
    ---------
    epochedt= copy of array obtained form mne epochs . epochs['stim_standard']) (if epochs from different events otherwise  np.array(epochs) + epoched.copy())
    returns array of epoched data with every other one inverted
    """

    for i in range(len(epochedt)):
        if not i % 2:
            epochedt[i, :, :] = -epochedt[i, :, :]

    return epochedt


def filtering(raw, notch=None, highpass=None, lowpass=None,
              fir_window="hamming", fir_design="firwin"):
    """
    base filter funciton

    Parameters
    ----------
    raw : raw mne object
    notch : list, single freqs to notch filter,
    other args self explanatory

    Returns
    -------
    raw : raw object filtered

    """

    _out_folder = 'C:\\Users\PC\Desktop\mea_anal'

    #   fig, ax = plt.subplots(2, sharex=True, sharey=True)
    #   fig.suptitle("Power Spectral Density")
    #   ax[0].set_title("before removing power line noise")
    #  ax[1].set_title("after removing power line noise")
    #   ax[1].set(xlabel="Frequency (Hz)", ylabel="μV²/Hz (dB)")
    #   ax[0].set(xlabel="Frequency (Hz)", ylabel="μV²/Hz (dB)")
    #  raw.plot_psd(average=True, area_mode=None, ax=ax[0], show=False)
    if notch is not None:  # notch filter at 50 Hz and harmonics
        raw.notch_filter(freqs=notch, fir_window=fir_window,
                         fir_design=fir_design)
    if lowpass is not None:  # lowpass filter at 50 Hz
        raw.filter(h_freq=lowpass, l_freq=None, fir_window=fir_window,
                   fir_design=fir_design)
    if highpass is not None:  # lowpass filter at 50 Hz
        raw.filter(h_freq=None, l_freq=highpass, fir_window=fir_window,
                   fir_design=fir_design)
    #   raw.plot_psd(average=True, area_mode=None, ax=ax[1], show=False)
    #   fig.tight_layout()
    #  fig.savefig(_out_folder/Path("remove_power_line_noise.pdf"), dpi=800)

    return raw

def dprimer(epochsdatasinglechannel, noiseref=[0,200]):
    "this function takes a 2d array of epochs, calculates average and std of the signal in a pre-stimulus epoch, and uses it to calculate dprime against a slider (avg and std ofdatapoint across trials)"
    dprime=np.zeros((len(epochsdata[0])))
    noiseavg=epochsdatasinglechannel[:,noiseref[0]:noiseref[1]].mean()
    noisestd=epochsdatasinglechannel[:,noiseref[0]:noiseref[1]].std()
    for i in range(0,len(epochsdatasinglechannel[0])):
        dprime[i]=abs(epochsdatasinglechannel[:,i].mean()-noiseavg)/np.sqrt((epochsdatasinglechannel[:,i].std()+noisestd)/2)
    return dprime
def bootstrapcore(epochsdatasinglechannel, function,functionarguments,trialstep=50,extract=1000):
        """"
        this function takes a 2d array of epochs and cycles through the trials taking a trial pool of the size of trialstep, increasing by trialstep
        (e.g. trialstep =50 the trial pool will be the first 50,100,150,200, 250 ecc trials successively)
        each pool is then sampled, for the same amount of trials, an amount of times equal to the extract input.
        so from the 50 trials pool, with extract= 1000, 1000 arrays of 50 trials, extracted from trialpool with replacement
        these arrays are each fed to an input function, the output of which is stored in a list
        """""
        storage=list()
        samplespool = epochsdatasinglechannel
        for r in range(0, int((len(epochsdatasinglechannel) / trialstep))):  # trial loop
            samplespoolt=samplespool[:(r+1) * trialstep,:]

            for bb in range(0, extract):
                inextract = np.random.randint(len(samplespoolt), size=(r+1) * trialstep)

                storage.append(function(samplespoolt[inextract],functionarguments))

        return storage



def Bootsrapper_snr(epochsdata,startrms=260,endrms=510,startnoiserms= 0, endnoiserms=250, startpeak=260, endpeak= 333 ,negpeak=333, negendpeak=420,secondpeak=390, secondendpeak=520, negsecondpeak=560, negsecondendpeak=810, noisepeak=0, noisendpeak=220 ):
    """
    this function calculates bootstrap snr in 50 trials increments. Also returns average trace

    Parameters
    ----------
    epochsdata : epoch data from mne epochs._data

    Returns
    -------
    all are per trial per channel
    sem : array, bootstrapped standard error of the mean snr (CHECK THIS)
    std : array, bootstrapped standard dev of snr.
    bpmean : ARRAY
        bootstrapped mean snr.
    confidencel :array
        confidence low snr.
    confidenceh : array 
        confidence high snr.

    """
    import scipy as sp
    from scipy import stats
    extract = 1000
    # snr calculated with rms, with sem and std
    bpsemrmssnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpstdrmssnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bprmssnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    # snr calculated with peak, p1 to n1, with sem and std
    bpsempeaksnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpstdpeaksnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bppeaksnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    # snr calculated with peak, n1 to p2, with sem and std
    bpsempeaktwosnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpstdpeaktwosnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bppeaktwosnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
# confidence intervals for all the above, at 5%
    bpconfidencelrmssnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpconfidencehrmssnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpconfidencelpeaksnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpconfidencehpeaksnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpconfidencelpeaktwosnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpconfidencehpeaktwosnr = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
# actual rms and peak values (not snr) bootstrapped
    bpsemrms = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpstdrms = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bprms = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpsempeakamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpstdpeakamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bppeakamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpsempeaktwoamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpstdpeaktwoamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bppeaktwoamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpconfidencelrms = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpconfidencehrms = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpconfidencelpeakamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpconfidencehpeakamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpconfidencelpeaktwoamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpconfidencehpeaktwoamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    # actual rms and peak values (not snr, not bootstrapped)
    rms = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    peakamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    peaktwoamp = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    #peak noise
    bpsempeaknoise= np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpstdpeaknoise= np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bppeaknoise = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))

    bpconfidencelpeaknoise = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))
    bpconfidencehpeaknoise = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1])))


    # we also want to store the actual trace with sem std and conf int

    # real trace, with sem and std
    trace = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1]), len(epochsdata[0,0,:])))
    stdtrace = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1]), len(epochsdata[0,0,:])))
    semtrace = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1]), len(epochsdata[0,0,:])))

    # confidence intervals for all the above, at 5%


    confidenceltrace = np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1]), len(epochsdata[0,0,:])))
    confidencehtrace= np.zeros((int(len(epochsdata[:, 0]) / 50), len(epochsdata[1]), len(epochsdata[0,0,:])))
    for u in range(0, len(epochsdata[1])):  # channel loop
        print("----------------------------", u)#track what channel iteration you are at
        samplespool = epochsdata[:, u, :]
        for r in range(0, int((len(epochsdata) / 50))):  # trial loop
            #here we save the real trace

            trace[r,u,:]=epochsdata[:(r+1)*50, u, :].mean(0)
            semtrace[r, u,:] = sp.stats.sem(epochsdata[:(r+1)*50, u, :], 0)
            stdtrace[r, u,:] = np.std(epochsdata[:(r+1)*50, u, :], 0)
            sortedtrace = np.sort(epochsdata[:(r+1)*50, u, :], 0)
            confidenceltrace[r, u,:] = sortedtrace[int((r+1)*50*5/100),:]
            confidencehtrace[r, u,:] = sortedtrace[int((r+1)*50*95/100),:]

            #and the real values of peaks and rms
            rms [r,u]= np.sqrt(np.mean(epochsdata[:(r+1)*50, u, :].mean(0)[startrms:endrms] ** 2))

            peakamp [r,u],gg= peak_minmax(epochsdata[:(r+1)*50, u, :].mean(0), pz=[startpeak, endpeak], nu=[negpeak, negendpeak])

            peaktwoamp [r,u],gg  = peak_minmax(epochsdata[:(r+1)*50, u, :].mean(0), pz=[secondpeak, secondendpeak], nu=[negsecondpeak, negsecondendpeak])



            #real trace done

            print(r) #track which channel amount iteration you are at

            #empty arrays for extraction of rms, peakamp 1&2 and repsective snrs
            extractedrmssnr = np.zeros(
                (extract))
            extractedpeaksnr = np.zeros(
                (extract))
            extractedpeaktwosnr = np.zeros(
                (extract))
            extractedpeakamp = np.zeros(
                (extract))
            extractedpeaknoise = np.zeros(
                (extract))
            extractedpeaktwoamp = np.zeros(
                (extract))
            extractedrms = np.zeros(
                (extract))
            samplespoolt=samplespool[:(r+1) * 50,:]
            for bb in range(0, extract):
                inextract = np.random.randint(len(samplespoolt), size=(r+1) * 50)
                #now we get our peak values and rms values
                extractedpeakamp[bb],  c= peak_minmax(samplespool[inextract].mean(0), pz=[startpeak, endpeak], nu=[negpeak, negendpeak]) #second output islatency of peak
                extractedpeaktwoamp[bb], cl= peak_minmax(samplespool[inextract].mean(0), pz=[secondpeak, secondendpeak], nu=[negsecondpeak, negsecondendpeak])
                extractedpeaknoise[bb],  cn = peak_minmax(samplespool[inextract].mean(0), pz=[noisepeak, noisendpeak], nu=[noisepeak, noisendpeak]) #this we don t save as it is useless outside of snr and can be calculated from snr
                extractedrms[bb]= np.sqrt(np.mean(samplespool[inextract].mean(0)[startrms:endrms] ** 2))
                extractedrmssnr[bb] = extractedrms[bb] / np.sqrt(np.mean(
                    samplespool[inextract].mean(0)[startnoiserms:endnoiserms] ** 2))  # we want our pool of samples for extraction. the pool of samples is the total of the trials for that channel, so it goes outside the trial loop
                extractedpeaksnr[bb] = extractedpeakamp[bb] / extractedpeaknoise[bb]
                extractedpeaktwosnr[bb] = extractedpeaktwoamp[bb] / extractedpeaknoise[bb]



            # at this point we have 10k snr values in 6 extracted variables: we get values for their stats, that s per tial and per channel
            #for rms snr
            bpsemrmssnr[r, u] = sp.stats.sem(extractedrmssnr, 0)
            bpstdrmssnr[r, u] = np.std(extractedrmssnr, 0)
            bprmssnr  [r, u] = np.mean(extractedrmssnr, 0)
            
            sortedr = np.sort(extractedrmssnr, 0)
            bpconfidencelrmssnr [r, u] = sortedr[int(extract*5/100)]
            bpconfidencehrmssnr[r, u] = sortedr[int(extract*95/100)]
            
            #now for peak 1 snr
            bpsempeaksnr[r, u] = sp.stats.sem(extractedpeaksnr, 0)
            bpstdpeaksnr[r, u] = np.std(extractedpeaksnr, 0)
            bppeaksnr[r, u] = np.mean(extractedpeaksnr, 0)

            sortedr = np.sort(extractedpeaksnr, 0)
            bpconfidencelpeaksnr[r, u] = sortedr[int(extract * 5 / 100)]
            bpconfidencehpeaksnr[r, u] = sortedr[int(extract * 95 / 100)]
            
            #for peak two snr
            bpsempeaktwosnr[r, u] = sp.stats.sem(extractedpeaktwosnr, 0)
            bpstdpeaktwosnr[r, u] = np.std(extractedpeaktwosnr, 0)
            bppeaktwosnr[r, u] = np.mean(extractedpeaktwosnr, 0)

            sortedr = np.sort(extractedpeaktwosnr, 0)
            bpconfidencelpeaktwosnr[r, u] = sortedr[int(extract * 5 / 100)]
            bpconfidencehpeaktwosnr[r, u] = sortedr[int(extract * 95 / 100)]
            
            
            #now all of the above for the values alone (not snr)

            # for rms 
            bpsemrms[r, u] = sp.stats.sem(extractedrms, 0)
            bpstdrms[r, u] = np.std(extractedrms, 0)
            bprms[r, u] = np.mean(extractedrms, 0)

            sortedr = np.sort(extractedrms, 0)
            bpconfidencelrms[r, u] = sortedr[int(extract * 5 / 100)]
            bpconfidencehrms[r, u] = sortedr[int(extract * 95 / 100)]

            # now for peak 1 amp
            bpsempeakamp[r, u] = sp.stats.sem(extractedpeakamp, 0)
            bpstdpeakamp[r, u] = np.std(extractedpeakamp, 0)
            bppeakamp[r, u] = np.mean(extractedpeakamp, 0)

            sortedr = np.sort(extractedpeakamp, 0)
            bpconfidencelpeakamp[r, u] = sortedr[int(extract * 5 / 100)]
            bpconfidencehpeakamp[r, u] = sortedr[int(extract * 95 / 100)]

            # for peak two amp
            bpsempeaktwoamp[r, u] = sp.stats.sem(extractedpeaktwoamp, 0)
            bpstdpeaktwoamp[r, u] = np.std(extractedpeaktwoamp, 0)
            bppeaktwoamp[r, u] = np.mean(extractedpeaktwoamp, 0)

            sortedr = np.sort(extractedpeaktwoamp, 0)
            bpconfidencelpeaktwoamp[r, u] = sortedr[int(extract * 5 / 100)]
            bpconfidencehpeaktwoamp[r, u] = sortedr[int(extract * 95 / 100)]
            
            # now for noise amp
            bpsempeaknoise[r, u] = sp.stats.sem(extractedpeaknoise, 0)
            bpstdpeaknoise[r, u] = np.std(extractedpeaknoise, 0)
            bppeaknoise[r, u] = np.mean(extractedpeaknoise, 0)

            sortedr = np.sort(extractedpeaknoise, 0)
            bpconfidencelpeaknoise[r, u] = sortedr[int(extract * 5 / 100)]
            bpconfidencehpeaknoise[r, u] = sortedr[int(extract * 95 / 100)]


    return     bpsemrmssnr ,bpstdrmssnr,bprmssnr,bpsempeaksnr,bpstdpeaksnr,bppeaksnr,bpsempeaktwosnr,bpstdpeaktwosnr,bppeaktwosnr,bpconfidencelrmssnr,bpconfidencehrmssnr,bpconfidencelpeaksnr,bpconfidencehpeaksnr,bpconfidencelpeaktwosnr,bpconfidencehpeaktwosnr,bpsemrms,bpstdrms,bprms,bpsempeakamp,bpstdpeakamp,bppeakamp,bpsempeaktwoamp,bpstdpeaktwoamp,bppeaktwoamp,bpconfidencelrms,bpconfidencehrms,bpconfidencelpeakamp,bpconfidencehpeakamp,bpconfidencelpeaktwoamp,bpconfidencehpeaktwoamp,rms,peakamp,peaktwoamp,trace,stdtrace,semtrace


def run_pipeline_n1(folder="Z:\\Alessandro_Braga\\MEA data february\\feb_n1_2021-02-26_16-35-13_m0002\\Record Node 101",session=1,lowpass=100, notch=50, reject_criteria = dict(eeg=300e-6), tmin=-0.04 ):
    import time
    if reject_criteria== "auto":
        reject_criteria=None
        autoreject= 1
    import mne  # this function cycles through all recordings for one recording session. b input is to tell where recordings'number starts from
    import pickle
    from glob import glob
    Folder = folder
    Files = sorted(glob(Folder + '/**/*.dat', recursive=True))

    for Fi, File in enumerate(Files):
        recording = Fi + 1  # 2*second animal
        print('now processing recording' + str(Fi + 1))

        indexo, filenamex, soa, stimdur, mean, data, triggers, timestamps, sfreq, oldfreq = loadOep(folder, session=session,
                                                                                                    recording=recording,
                                                                                                    newfreq=5000)
        import matplotlib.pyplot as plt

        raw, events = MNEify(indexo, data, timestamps, triggers, sfreq, oldfreq=5000)
        #hack for february data to align
        if recording==1:
            events[:,0]+=165
        elif recording == 2:
            events[:, 0] += 185
        elif recording == 3:
            events[:, 0] += 45

        raw = filtering(raw, notch=notch, highpass=None, lowpass=None,
                        fir_window="hamming", fir_design="firwin")
        raw = filtering(raw, notch=None, highpass=None, lowpass=lowpass,
                        fir_window="hamming", fir_design="firwin")
        epochs = mne.Epochs(raw, events, tmin=tmin, reject=reject_criteria, baseline=None, preload=True, flat=None, proj=False,
                            reject_by_annotation=False)
        if autoreject:
            from autoreject import AutoReject  # import the autoreject module . all this stuff is functionally useless
            ar = AutoReject(n_interpolate=[3, 6, 12], random_state=42)
            epochs, reject_log = ar.fit_transform(epochs, return_log=True)
        # ica = mne.preprocessing.ICA(n_components=0.99, method="fastica")
        # ica.fit(epochs_ar)
        # ica_sources = ica.get_sources(epochs_ar)
        # rr = ica_sources._data
        # rrr = rr[:, :, :].mean(0)
        # byo = np.zeros((len(rrr)))
        # byor = np.zeros((len(rrr)))
        # for i in range(0, len(rrr)):
        #     byo[i] = (max(rrr[i, 250:500]) - min(rrr[i, 250:500]))
        #     byor[i] = (max(rrr[i, 500:]) - min(rrr[i, 500:]))
        # epochs_ica = ica.apply(epochs_ar, exclude=np.argwhere(byo == min(byo))[0])

        epochsdata = epochs._data
        bpsemrmssnr, bpstdrmssnr, bprmssnr, bpsempeaksnr, bpstdpeaksnr, bppeaksnr, bpsempeaktwosnr, bpstdpeaktwosnr, bppeaktwosnr, bpconfidencelrmssnr, bpconfidencehrmssnr, bpconfidencelpeaksnr, bpconfidencehpeaksnr, bpconfidencelpeaktwosnr, bpconfidencehpeaktwosnr, bpsemrms, bpstdrms, bprms, bpsempeakamp, bpstdpeakamp, bppeakamp, bpsempeaktwoamp, bpstdpeaktwoamp, bppeaktwoamp, bpconfidencelrms, bpconfidencehrms, bpconfidencelpeakamp, bpconfidencehpeakamp, bpconfidencelpeaktwoamp, bpconfidencehpeaktwoamp, rms, peakamp, peaktwoamp, trace, stdtrace, semtrace=Bootsrapper_snr(epochsdata)
        savedict={"bpsemrmssnr":bpsemrmssnr,"bpstdrmssnr":bpstdrmssnr,"bprmssnr":bprmssnr,"bpsempeaksnr":bpsempeaksnr,"bpstdpeaksnr":bpstdpeaksnr,"bppeaksnr":bppeaksnr,"bpsempeaktwosnr":bpsempeaktwosnr,"bpstdpeaktwosnr":bpstdpeaktwosnr,"bppeaktwosnr":bppeaktwosnr,"bpconfidencelrmssnr":bpconfidencelrmssnr, "bpconfidencehrmssnr":bpconfidencehrmssnr, "bpconfidencelpeaksnr":bpconfidencelpeaksnr,"bpconfidencehpeaksnr":bpconfidencehpeaksnr,"bpconfidencelpeaktwosnr":bpconfidencelpeaktwosnr,"bpconfidencehpeaktwosnr":bpconfidencehpeaktwosnr, "bpsemrms":bpsemrms,"bpstdrms":bpstdrms,"bprms":bprms,"bpsempeakamp":bpsempeakamp, "bpstdpeakamp":bpstdpeakamp, "bppeakamp":bppeakamp, "bpsempeaktwoamp":bpsempeaktwoamp, "bpstdpeaktwoamp":bpstdpeaktwoamp, "bppeaktwoamp":bppeaktwoamp, "bpconfidencelrms":bpconfidencelrms, "bpconfidencehrms":bpconfidencehrms, "bpconfidencelpeakamp":bpconfidencelpeakamp, "bpconfidencehpeakamp":bpconfidencehpeakamp, "bpconfidencelpeaktwoamp":bpconfidencelpeaktwoamp, "bpconfidencehpeaktwoamp":bpconfidencehpeaktwoamp, "rms":rms, "peakamp":peakamp, "peaktwoamp":peaktwoamp, "trace":trace,"stdtrace":stdtrace,"semtrace":semtrace}
        with open(folder[45:70]+'_bootstrapOGsnr_&_amp_recording'+str(recording)+'_'+str(time.time())[:8]+'.pickle', 'wb') as handle:
            pickle.dump(savedict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(folder[45:70]+'_reject_log_recording'+str(recording)+'_'+str(time.time())[:8]+'.pickle', 'wb') as handle:
            pickle.dump(reject_log, handle, protocol=pickle.HIGHEST_PROTOCOL)
        # indexo, filenamex, soa, stimdur, mean, data, datar, triggers, timestamps, sfreq = loadOep(folder, session,
        #                                                                                           recording)  # get data
        #
        # if len(triggers) < 100:  # (to avoid processing recordings that are botched)
        #     continue
        # raw, events = MNEify(indexo, data, timestamps, triggers, sfreq)
        # raw = filtering(raw, notch=50, highpass=None, lowpass=None,
        #                 fir_window="hamming", fir_design="firwin")
        # raw = filtering(raw, notch=None, highpass=None, lowpass=100,
        #                 fir_window="hamming", fir_design="firwin")
        # reject = dict(
        #     eeg=280e-6,  # (EEG channels)
        #
        # )
        # if np.array(
        #         indexo).any():  # hereby start processing experiments with index (experiments in which trials have different conditions/stimuli)
        #
        #     event_dict = {'stim_standard': 1, 'omission': 0}
        #     epochs = mne.Epochs(raw, events, event_id=event_dict, tmin=-0.06, reject=reject, baseline=None,
        #                         preload=True, flat=None, proj=False, reject_by_annotation=False)
        #     while len(epochs.events) < len(
        #             events) * 0.8:  # this while loop is used to reject bad channels, that is channels that force dropping toomany epochs.
        #
        #         bd = epochs.ch_names
        #         prol = np.sort(data, 0)
        #         ro = np.zeros((4000))
        #
        #         for dd in range(1, 4000):
        #             rr = np.zeros((len(data[0])))
        #             for i in range(0, len(data[0])):
        #                 rr[i] = np.mean(prol[:dd, i] - prol[-dd:, i])
        #             ro[dd - 1] = (np.argwhere(rr == min(rr)))
        #         for u, uu in enumerate(bd):
        #             rr[u] = sum(ro == u)
        #         todrop = np.argwhere(rr == max(rr))
        #         a = [str(bd[todrop[0, 0]])]
        #         raw.info['bads'] = a
        #         print(str(bd[todrop[0, 0]]))
        #         raw.pick_types(eeg=True, exclude='bads')
        #
        #         epochs = mne.Epochs(raw, events, tmin=-0.06, event_id=event_dict, reject=reject, baseline=None,
        #                             preload=True, flat=None, proj=False, reject_by_annotation=False, verbose='WARNING')
        #         data = np.delete(data, todrop[0, 0], 1)
        #     epochedstim = np.array(epochs['stim_standard'])
        #     epochedomis = np.array(epochs['omission'])
        #
        # else:  # ereby processing of experiments in which laltrials are the same(no index)
        #
        #     epochs = mne.Epochs(raw, events, tmin=-0.06, reject=reject, baseline=None, preload=True, flat=None,
        #                         proj=False, reject_by_annotation=False)
        #     # dropping=[]
        #     while len(epochs.events) < len(
        #             events) * 0.8:  # this loop to identify bad channels and avoid dropping too many epochs.
        #         bd = epochs.ch_names
        #
        #         prol = np.sort(data, 0)
        #         ro = np.zeros((4000))
        #
        #         for dd in range(1, 4000):
        #             rr = np.zeros((len(data[0])))
        #             for i in range(0, len(data[0])):
        #                 rr[i] = np.mean(prol[:dd, i] - prol[-dd:, i])
        #             ro[dd - 1] = (np.argwhere(rr == min(rr)))
        #         for u, uu in enumerate(bd):
        #             rr[u] = sum(ro == u)
        #         todrop = np.argwhere(rr == max(rr))
        #         a = [str(bd[todrop[0, 0]])]
        #         # dropping.append(a[0])
        #         raw.info['bads'] = a
        #         print(str(bd[todrop[0, 0]]))
        #         raw.pick_types(eeg=True, exclude='bads')
        #         epochs = mne.Epochs(raw, events, tmin=-0.06, reject=reject, baseline=None, preload=True, flat=None,
        #                             proj=False, reject_by_annotation=False, verbose='WARNING')
        #         data = np.delete(data, todrop[0, 0], 1)
        #     epoched = np.array(epochs)
        #     evoked = epochs.average()


#
# run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\SOAMMNSOA_2020-12-09_15-43-32_m0004\Record Node 101",session=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=1)#
# run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\SOAMMNSOA_2020-12-01_14-09-21_m0003\Record Node 101",session=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=2)
# run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\SOANONSOA_2020-11-16_15-14-52_m0002\Record Node 101",session=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=1)
#run_pipeline_n1(folder = "Z:\\Alessandro_Braga\\MEA data february\\feb_n1_2021-02-26_16-35-13_m0002\Record Node 101",session=1)#
# folder="Z:\\Alessandro_Braga\\MEA data february\\feb_n1_2021-02-26_16-35-13_m0002\Record Node 101"
# indexo,filenamex,soa, stimdur,mean,data,datar,triggers,timestamps,sfreq,oldfreq=loadOep(folder, session=1, recording=3)#get data
#
# folder = "Z:\\Alessandro_Braga\\MEA data february\\feb_n1_2021-02-26_16-35-13_m0002\Record Node 101"
# indexo, filenamex, soa, stimdur, mean, data, triggers, timestamps, sfreq, oldfreq = loadOep(folder, session=1,
#                                                                                             recording=2, newfreq=5000)  # get data
# import matplotlib.pyplot as plt
#
# raw, events = MNEify(indexo, data, timestamps, triggers, sfreq, oldfreq=5000)
# raw = filtering(raw, notch=50, highpass=None, lowpass=None,
#                 fir_wi.99, method="fastica")
# ica.fit(epochs_ar)
# ica_sources = ica.get_sources(epochs_ar)ndow="hamming", fir_design="firwin")
# raw = filtering(raw, notch=None, highpass=None, lowpass=100,
#                 fir_window="hamming", fir_design="firwin")
# epochs = mne.Epochs(raw, events, tmin=-0.04, reject=None, baseline=None, preload=True, flat=None, proj=False,
#                     reject_by_annotation=False)
# from autoreject import AutoReject # import the autoreject module .
# ar=AutoReject(n_interpolate=[3,6,12], random_state=42)
# epochs_ar, reject_log = ar.fit_transform(epochs, return_log=True)
# ica = mne.preprocessing.ICA(n_components=0
#
#
# rr=ica_sources._data
# rrr=rr[:,:,:].mean(0)
# byo = np.zeros((len(rrr)))
# for i in range(0, len(rrr)):
#     byo[i] = (max(rrr[i, 250:500]) - min(rrr[i, 250:500]))
# epochs_ica = ica.apply(epochs_ar, exclude=np.argwhere(byo==min(byo))[0])
#
# epochsdata = epochs_ica._data

# #
# # # dropping=[]
# #
# # epochsdata = epochs._data
# evoked=epochs.average()
# evoked.plot()
# ts_args = dict(time_unit='s')
# topomap_args = dict(sensors=True, time_unit='s')
# evoked.plot_joint(title='right auditory',
#                         ts_args=ts_args, topomap_args=topomap_args)

"""
all these functions rely on filenames and folders retaining the structure set by oeps. h5 files of a certain recording-experiment must be put in the corresponding recordingN (e.g. recording5) folder.
"""

def loadOep(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\SOAMMNSOA_2020-12-01_14-09-21_m0003\Record Node 101", session=1,newfreq=500, recording=10, mainevent=8):
    """

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
    omiss=False
    indexo=[]#it stays [] if all trials are identical and there is no index
    front_left=['B6' ,'B5' ,'C6' ,'C5']#for future use
    channelmap=np.array([[16 ,18 ,15 ,17],[14 ,12 ,13 ,11],[10,0,9,7],[20,29,19,22],[6,3,5,2],[23,27,24,26],[8,1,21,28],[28,1,4,25]])#see below
    
    front_right=['D6' ,'D5' ,'E6' ,'E5']
    front_righti=[14 ,12 ,13 ,11]
    media_right=['D4' ,'F3' ,'E4' ,'E3']
    media_righti=[10,0,9,7]
    media_left=['B4' ,'A3' ,'C4' ,'B3']
    media_lefti=[20,29,19,22]
    poste_right=['F1' ,'F2' ,'E1' ,'E2']
    poste_righti=[6,3,5,2]
    poste_left=['A1' ,'B2' ,'B1' ,'A2']
    poste_lefti=[23,27,24,26]
    cent=['D3' ,'C2' ,'C3' ,'D2']
    centi=[8,1,21,28]
    poste_cent=['D1' ,'C2' ,'C1' ,'D2']
    poste_centi=[28,1,4,25]
    _out_folder = 'C:\\Users\PC\Desktop\mea_anal'

    Folder=folder #folder where the recordings are
    session=session#usually 1. it is the overall recording session
    Recording=recording# various recordings per session
    Processor=None
    Files = sorted(glob(Folder+'/**/*.dat', recursive=True))#get all data files in folder
    InfoFiles = sorted(glob(Folder+'/*/*/structure.oebin'))


    Data, Rate = {}, {}
    for F,File in enumerate(Files):#go through, select the recording  from input.
        
        
        Exp, Rec, _, Proc = File.split('\\')[-5:-1]
        Exp = str(int(Exp[10:])-1)
        Rec = str(int(Rec[9:])-1)
        Proc = Proc.split('.')[0].split('-')[-1]
        if '_' in Proc: Proc = Proc.split('_')[0]
        
        if Proc not in Data.keys(): Data[Proc], Rate[Proc] = {}, {}
        
        if session:
            if int(Exp) != session-1: continue
        
        if Recording:
            if int(Rec) != Recording-1: continue
        
        if Processor:
            if Proc != Processor: continue
        
        print('Loading recording', int(Rec)+1, '...')
        if Exp not in Data[Proc]: Data[Proc][Exp] = {}
        Data[Proc][Exp][Rec] = np.memmap(File, dtype='int16', mode='c')
        
        
        Info = literal_eval(open(InfoFiles[F]).read())
        ProcIndex = [Info['continuous'].index(_) for _ in Info['continuous']
                     if str(_['recorded_processor_id']) == '101'][0]
        
        ChNo = Info['continuous'][ProcIndex]['num_channels']
        if Data[Proc][Exp][Rec].shape[0]%ChNo:
            print('Rec', Rec, 'is broken')
            del(Data[Proc][Exp][Rec])
            continue
        
        SamplesPerCh = Data[Proc][Exp][Rec].shape[0]//ChNo
        Data[Proc][Exp][Rec] = Data[Proc][Exp][Rec].reshape((SamplesPerCh, ChNo))
        Rate[Proc][Exp] = Info['continuous'][ProcIndex]['sample_rate']
        f=Rate[Proc][Exp]
        oldfreq=f
        FolderEv=Folder+'\\experiment'+str(session)+'\\recording'+str(recording)+'\\events\Rhythm_FPGA-100.0\\TTL_1'#get folder with events
        FilesEv = sorted(glob(FolderEv+'\\*.npy', recursive=True))#get events
        ttll=np.load(FilesEv[0])
            
        ttl=np.argwhere(ttll==mainevent)#get trial event.
        tstp=np.load(FilesEv[3])#load timestamps
        data=Data[Proc][Exp][str(recording-1)]
        chunk=int((len(Data[Proc][Exp][str(recording-1)])/len(ttl))/(int(f/newfreq) ))#size of trial in datapoints


        filenamex=File#to set aside the filename of the actual recording we are working with after we exit the loop through all FIles for that session

        u=[]
   
   




    Data=[]

    
    myfile=filenamex[:-3]+str(newfreq)+'_dwnsmpl.npy'#this is to average together channels as indicated in the top list (channelmap). looked weird when i actually used these data TODO look into averaging together neighboring channels to go from 30 to 8
  
    myfilec=filenamex[:-3]+'_dwnsmpl_all_channels.npy'#
    from pathlib import Path

    my_file = Path(myfile)
    if my_file.is_file():
        print('load downsampled file')

        data=np.load(myfilec)
        datar=np.load(myfile)
        f=newfreq
    else:
        print('saving downsampled file')#if frequency is not newfreq hz data are downsampled to newfreq
        if f!=newfreq:
            data= data[::int(f/newfreq),:32]#downsample from 30k
            f=newfreq
            
        data=np.delete(data,[7,24],1)#get rid of non-channels
        datar=np.zeros((len(data),len(channelmap)))
        for i,d in enumerate(channelmap):
            datar[:,i]=data[:,d].mean(1)
            
        data=(data*0.195/1000000)#convert bits to microvolt and microvolt to volt. "The raw data are saved as signed 16-bit integers, in the range -32768 to 32767. They don’t have a unit. To convert to microvolts, just  multiply by 0.195, as you noted. This scales the data into the range ±6.390 mV, with 0.195 µV resolution (Intan chips have a ±5 mV input range).""
        datar=(datar*0.195/1000000)#convert bits to microvolt and microvolt to volt. THIS GETS DATA IN ACTUAL MEMORY; TAKES TIME
        data=data-data.mean(0)#baseline
        datar=datar-datar.mean(0)
        np.save(myfile, datar)
        np.save(myfilec, data)
    
    h5file=sorted(glob(Folder+'/experiment'+str(session)+'/recording'+str(recording)+'/**/*.h5', recursive=True))#this is the file from the tdt side. they are put by hand in the recordingn folder after data collection
    import tables
    try:
        ff=tables.open_file(h5file[0])
        soa=((ff.root._v_attrs.Exp_Settings_SOA/ff.root._v_attrs.Exp_Settings_sampling_freq))#SOA for filenames. rounding issue here 
        stimdur=((ff.root._v_attrs.Exp_Settings_stim_dur/ff.root._v_attrs.Exp_Settings_sampling_freq))#stim duration for filenames
        if h5file[0][-14]=='N':# checks if the recording is a MMN one, in which case gets the index.
            omiss=True  
            indexo=Omitter(h5file[0])#get index of omissions from recorded audio signal. this allows to check for discrepancies from index (there should be none if the computer is not used during recording)
    except Exception:#in case the h5 file is damaged or does not exist.
        soa=((len(datar)/(newfreq**2)))
        stimdur=0.109 
    index=tstp[ttl]/(oldfreq/newfreq)#downsample timestamps
    index=index-index[0]
    if index[-1]+chunk>len(datar):#hack  if recording happens to be a datapoint too short last trial is deleted.
            index=np.delete(index,len(index)-1)    
            deletable=np.argwhere((tstp-tstp[0])/int(oldfreq/newfreq)+chunk>len(datar))
            tstp=np.delete(tstp,deletable)
            ttll=np.delete(ttll,deletable)
    u=[]  #these following loops are to average epochs , the purpose being having quick way to look at erp after recording session
    for iz in range (0,len(datar[1,:])):
        m=[]
     
        for i in range(0,len(index)):
           m.append(datar[int(index[i]):int(index[i])+chunk,iz])
    
        u.append(sum(m)/len(m)-sum(sum(m))/(len(m)*len(sum(m))))
        
        v=[]
    for iz in range (0,len(data[1,:])):
        n=[]
        

        for i in range(0,len(index)):
           n.append(data[int(index[i]):int(index[i])+chunk,iz])

        v.append(sum(n)/len(n)-sum(sum(n))/(len(n)*len(sum(n))))


    mean=v#this is the epoched and averaged data
    triggers=ttll
    sfreq=f
    timestamps=tstp
    return indexo,filenamex,soa,stimdur,mean,data,datar,triggers,timestamps,sfreq,oldfreq
#indexo: stimulus index if mmn from h5. filename. stimulus onset async from h5. stimulus duration from h5. mean=averaged epoched traces, data = downsampled data, datar= downsampled data iwth channels averaged in groups (left front ecc). triggers= trigger events index for timestamps in teimestams, sfrequ is sampling freq should always be newfreq at this pint
    

def Omitter(h5file):
    """
    this method goes through the audio recording and checks which trials are silent. this is a temporary solution for omission MMN as in two recordings the index was not followed (intensity switching logic is time dependent and can fail. Should be ok now)
    """
    omiss:True      
    
    filename=h5file
    
    
    import tables

    import numpy as np
    import numpy


    
    ff=tables.open_file(filename)
    leftn1=[]
    for group in ff.walk_groups("/"):
         for array in ff.list_nodes(group, classname='Array'):
             leftn1.append(array[:,1])
    leftn1=np.array(leftn1)*50000
    leftarr=np.zeros((ff.root._v_attrs.Exp_Settings_reps,int(ff.root._v_attrs.Exp_Settings_SOA/10)))
    for i in range(0,ff.root._v_attrs.Exp_Settings_reps):
         leftarr[i]=numpy.array(leftn1)[0,(i*int(ff.root._v_attrs.Exp_Settings_SOA/10)) : (i+1)*int(ff.root._v_attrs.Exp_Settings_SOA/10)]
    dd=np.max(leftarr[:,375:625],1)
    uu=dd>100
    return uu
def scaler(my_diction):  
    total = 0  
     
    for j in my_diction:  
        my_diction[j] = (my_diction[j])/1000
    return my_diction
def MNEify(indexo,datar,timestamps,triggers,sfreq,oldfreq,trialevent=8): 
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

    Returns
    -------
    raw objects and event array to make epochs
    """
            
        import mne
        datad=np.transpose(datar[:,:len(datar[1,:])])#transpose raw data as per mne requirment
        timstp=timestamps-timestamps[0]#zero timestamps to start of recording
        eventstsp=timstp[np.argwhere(triggers==trialevent)]#select timestamps of each trial event
        ggg=np.zeros((len(eventstsp),1))
        dd=np.ones((len(eventstsp),1))
        if np.array(indexo).any(): #fix discrepancies between h5file index and timestamps
            dd=indexo.astype(int)
            dd=np.reshape(dd,(len(dd),1))
            if len(dd)>len(eventstsp):
                dd=dd[:len(eventstsp)]
            elif len(dd)<len(eventstsp):
                eventstsp=eventstsp[:len(dd)]
                ggg=ggg[:len(dd)]
                
#     tim=tim-tim[1843]#3000 datapoint at 48828.125hz recorded at 30000hz

# unfilt=np.delete(unfilt,[7,24],0)
# filtt=np.delete(filtt,[7,24],0)
        if len(datad)==30:#channel names for all 30
            import pandas as pd

            df = pd.read_csv('channelsH32Mabs.txt')
            ch_names = df.name.to_list()

            pos = df[['z','x','y']].values
            dig_ch_pos = dict(zip(ch_names,pos))
            scaler(dig_ch_pos)
            nasion= [-0.007, 0, -0.001]
            lpa= [0.000, -0.005, -0.001]
            rpa=[0.000, 0.005, -0.001]
            montage = mne.channels.make_dig_montage(ch_pos=dig_ch_pos,nasion=nasion,lpa=lpa,rpa=rpa)
                        
            ch_types = ['eeg', 'eeg', 'eeg', 'eeg','eeg', 'eeg', 'eeg','eeg', 'eeg', 'eeg', 'eeg','eeg', 'eeg', 'eeg', 'eeg','eeg', 'eeg', 'eeg', 'eeg','eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg','eeg', 'eeg', 'eeg', 'eeg']
        else: #channel names for pooled channels
            ch_names = ['front_left','front_right','cent_right','cent_left','post_right','post_left','central','post_cent']
            ch_types = ['eeg', 'eeg', 'eeg', 'eeg','eeg', 'eeg', 'eeg','eeg']
        
            
        

        events=np.concatenate((eventstsp/(oldfreq/sfreq)+30,ggg,dd),1)#30 is fixed from rcx delay. here we get events in shape required from mne (1st column timestamps, second zeros, third ones)
        events=events.astype(int)
        info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
        raw = mne.io.RawArray(datad, info)
        raw=raw.set_montage(montage)
        events=events.astype(int)
        return raw,events
        
        
from mne import set_eeg_reference, events_from_annotations
from mne.epochs import Epochs
from autoreject import AutoReject, Ransac
import numpy as np
from mne.preprocessing.ica import ICA, corrmap, read_ica
import os
from pathlib import Path
from matplotlib import pyplot as plt, patches



        
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
    peak_idx = peakutils.indexes(data-base, thres=thresh, min_dist=min_dist)
    if peak_idx.size==0:
        peak_idx=np.argwhere(data==max(data))
    return peak_idx, base
        
def peak_minmax(trace,pz=[30,45],nu=[45,60]):
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

    bb=find_peaks(trace[pz[0]:pz[1]], min_dist=1, thresh=0.3, degree=None)
    tt=find_peaks(0-trace[nu[0]:nu[1]], min_dist=1, thresh=0.0, degree=None)
    pospeakindex=pz[0]+np.argwhere(trace==max(trace[pz[0]+bb[0]]))
    negpeak=(max(0-trace[nu[0]+tt[0]]))
    pospeak=(max(trace[pz[0]+bb[0]]))
    negpeakindex=nu[0]+np.argwhere(0-trace==max(0-trace[nu[0]+tt[0]]))
    return pospeak,negpeak,pospeakindex,negpeakindex
def noise_rms(epochedt,chan,tr):
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
    evokedchannel=np.array(epochedt)[0:tr,chan,:].mean(0)
    rms = np.sqrt(np.mean(evokedchannel**2))
    
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

def snr_peak_by_trial(epoched):
    """
        
     calculate snr in a number of different manners
    

    Parameters
    ----------
    epoched : array
         obtained as np.array(epochs['stim_standard']) (if epochs from different events otherwise  np.array(epochs) )


    Returns
    -------
    all arrays
    rmsnp=rms of n1 complex
    snr_before_event=rmsnp divided by rms of epoch same size pre event, 
    snr_later=as above but after as far as poss from event, 
    snr_avg_inverted=rmsnp divided by the rms calculated on the inverted epochs(see noise inverter func),
    peaksnr= snr calculated on peak amplitudes
    peaksnrwide,snr calculated on peaks from a wide epoch
    peakamp, amplitude peak to peak
    posindex,index of positive peak after event
    negindex= index of neg peak after event
    """
    epochedt=epoched.copy()
    snr_before_event=np.zeros((len(epoched[1]),len(epoched)))
    rmsnp=np.zeros((len(epoched[1]),len(epoched)))
    snr_later=np.zeros((len(epoched[1]),len(epoched)))
    snr_avg_inverted=np.zeros((len(epoched[1]),len(epoched)))
    bbbb=np.zeros((len(epoched[1]),len(epoched)))
    mmmm=np.zeros((len(epoched[1]),len(epoched)))
    posindex=np.zeros((len(epoched[1]),len(epoched)))
    negindex=np.zeros((len(epoched[1]),len(epoched)))
    lllll=np.zeros((len(epoched[1]),len(epoched)))
    ddddd=np.zeros((len(epoched[1]),len(epoched)))
    bbbbb=np.zeros((len(epoched[1]),len(epoched)))
    maxmin=np.zeros((len(epoched[1]),len(epoched)))
    mmmmm=np.zeros((len(epoched[1]),len(epoched)))
    inverted=np.sqrt(np.mean(noise_inverter(epochedt).mean(0)**2))
    for gg in range(0,len(epoched[1])):
        print(gg)   
        for ii in range(1,len(epoched)):

            rmsnp[gg,ii]=np.sqrt(np.mean(np.mean(epoched[0:ii,gg,30:60],0)**2))#for epoched at -0.6 tith 30 dp from event, comparing with before event
            maxmin[gg,ii]=np.max((np.mean(epoched[0:ii,gg,150:],0)))-np.min((np.mean(epoched[0:ii,gg,150:],0)))

            snr_before_event[gg,ii]=np.sqrt(np.mean(np.mean(epoched[0:ii,gg,30:60],0)**2))/np.sqrt(np.mean(np.mean(epoched[0:ii,gg,0:30],0)**2))#for epoched at -0.6 tith 30 dp from event, comparing with before event
            snr_later[gg,ii]=np.sqrt(np.mean(np.mean(epoched[0:ii,gg,30:60],0)**2))/np.sqrt(np.mean(np.mean(epoched[0:ii,gg,-40:-10],0)**2))#for epoched at -0.6 tith 30 dp from event, comparing with latest from event
            snr_avg_inverted[gg,ii]=np.sqrt(np.mean(np.mean(epoched[0:ii,gg,30:60],0)**2))/inverted#for epoched at -0.6 tith 30 dp from event, comparing with Schimmel(1967)
            bbbb[gg,ii],mmmm[gg,ii],posindex[gg,ii],negindex[gg,ii]=peak_minmax((np.mean(epoched[0:ii,gg,:],0)))#this saves Voltage of p0 and n1 and their positions
            bbbbb[gg,ii],mmmmm[gg,ii],lllll[gg,ii],ddddd[gg,ii]=peak_minmax((np.mean(epoched[0:ii,gg,:],0)),[0,15],[15,30])#this gets a pos and neg peak from the pre stim interval
    peakamp=(bbbb-mmmm)
    peaknoise=np.sqrt((bbbbb-mmmmm)**2)
    peaksnr=peakamp/peaknoise
    peaksnr=peakamp/peaknoise
    peaksnrwide=peakamp/maxmin
    
    return rmsnp,snr_before_event,snr_later,snr_avg_inverted,peaksnr,peaksnrwide,peakamp,posindex,negindex




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



    
def Bootsrapper_snr(epochsdata):
    """
    this function calculates bootstrap snr. it takes ages (it bootstraps every amount of trials for every channel)

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
    extract=1000
    #we need to store ntrial snr values, ntrial std, ntrial
    sem=np.zeros((len(epochsdata[:,0]),len(epochsdata[1])))
    std=np.zeros((len(epochsdata[:,0]),len(epochsdata[1])))
    bpmean=np.zeros((len(epochsdata[:,0]),len(epochsdata[1])))

    confidencel=np.zeros((len(epochsdata[:,0]),len(epochsdata[1])))
    confidenceh=np.zeros((len(epochsdata[:,0]),len(epochsdata[1])))
    



    for u in range(0,len(epochsdata[1])):#channel loop
        print(u)
        samplespool=epochsdata[:,u,:]
        for r in range(1,len(epochsdata)):# trial loop
            
            extracted=np.zeros((extract))#first we make empty array for extraction. this array needs to hold 10000 rms/rms of mean
            for bb in range (0,extract):
                      inextract=np.random.randint( len( samplespool), size=r)
                      extracted[bb]= np.sqrt(np.mean(samplespool[inextract].mean(0)[30:60]**2))    /  np.sqrt(np.mean(samplespool[inextract].mean(0)[0:30]**2)   )    #we want our pool of samples for extraction. the pool of samples is the total of the trials for that channel, so it goes outside the trial loop
            
            #at this point we have 10k snr values: we get values for their stats, that s per tial and per channel 
            sem[r,u]=sp.stats.sem(extracted,0)
            std[r,u]=np.std(extracted,0)
            bpmean[r,u]=np.mean(extracted,0)

            sortedr=np.sort(extracted,0)
            confidencel[r,u]=sortedr[60]
            confidenceh[r,u]=sortedr[940]  


    return sem,std,  bpmean, confidencel,confidenceh

def run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\TESTAUD_2020-11-16_15-14-52_m0002_1sandothers\Record Node 101",session=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=1):
    import mne#this function cycles through all recordings for one recording session. b input is to tell where recordings'number starts from
    global _out_folder
    _out_folder = out_folder
    from glob import glob
    Folder=folder
    Files = sorted(glob(Folder+'/**/*.dat', recursive=True))

    for Fi,File in enumerate(Files):
        recording=Fi+b#2*second animal
        print('now processing recording'+str(Fi+b))
        indexo,filenamex,soa, stimdur,mean,data,datar,triggers,timestamps,sfreq=loadOep(folder, session, recording)#get data
        
        if len(triggers)<100:# (to avoid processing recordings that are botched)
            continue
        raw,events=MNEify(indexo,data,timestamps,triggers,sfreq)
        raw=filtering(raw, notch=50, highpass=None, lowpass=None,
              fir_window="hamming", fir_design="firwin")
        raw=filtering(raw, notch=None, highpass=None, lowpass=100,
              fir_window="hamming", fir_design="firwin")
        reject = dict(
              eeg=280e-6, #  (EEG channels)

              )
        if np.array(indexo).any():                 #hereby start processing experiments with index (experiments in which trials have different conditions/stimuli)

            event_dict = {'stim_standard': 1, 'omission': 0}
            epochs = mne.Epochs(raw, events, event_id=event_dict, tmin= -0.06, reject=reject, baseline=None, preload=True, flat = None, proj=False,reject_by_annotation=False)
            while len(epochs.events)<len(events)*0.8:                 #this while loop is used to reject bad channels, that is channels that force dropping toomany epochs.

                bd=epochs.ch_names
                prol=np.sort(data,0)
                ro=np.zeros((4000))
                
                for dd in range(1,4000):
                    rr=np.zeros((len(data[0])))
                    for i in range(0,len(data[0])):
                        rr[i]=np.mean(prol[:dd,i]-prol[-dd:,i])
                    ro[dd-1]=(np.argwhere(rr==min(rr)))
                for u,uu in enumerate(bd):
                    rr[u]=sum(ro==u)
                todrop=np.argwhere(rr==max(rr))
                a=[str( bd[todrop[0,0]])]
                raw.info['bads']=a
                print(str( bd[todrop[0,0]]))
                raw.pick_types(eeg=True, exclude='bads')
                
                epochs = mne.Epochs(raw, events, tmin= -0.06,event_id=event_dict, reject=reject, baseline=None, preload=True, flat = None, proj=False,reject_by_annotation=False, verbose='WARNING')
                data=np.delete(data,todrop[0,0],1)
            epochedstim=np.array(epochs['stim_standard']) 
            epochedomis=np.array(epochs['omission']) 

        else:                 #ereby processing of experiments in which laltrials are the same(no index)
           
            epochs = mne.Epochs(raw, events, tmin= -0.06, reject=reject, baseline=None, preload=True, flat = None, proj=False,reject_by_annotation=False)
            #dropping=[]
            while len(epochs.events)<len(events)*0.8:                #this loop to identify bad channels and avoid dropping too many epochs.
                bd=epochs.ch_names

            
                
            
                prol=np.sort(data,0)
                ro=np.zeros((4000))
                
                for dd in range(1,4000):
                    rr=np.zeros((len(data[0])))
                    for i in range(0,len(data[0])):
                        rr[i]=np.mean(prol[:dd,i]-prol[-dd:,i])
                    ro[dd-1]=(np.argwhere(rr==min(rr)))
                for u,uu in enumerate(bd):
                    rr[u]=sum(ro==u)
                todrop=np.argwhere(rr==max(rr))
                a=[str( bd[todrop[0,0]])]
                #dropping.append(a[0])
                raw.info['bads']=a
                print(str( bd[todrop[0,0]]))
                raw.pick_types(eeg=True, exclude='bads')
                epochs = mne.Epochs(raw, events, tmin= -0.06, reject=reject, baseline=None, preload=True, flat = None, proj=False,reject_by_annotation=False, verbose='WARNING')
                data=np.delete(data,todrop[0,0],1)
            epoched=np.array(epochs) 
            evoked=epochs.average()

#
#run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\SOAMMNSOA_2020-12-09_15-43-32_m0004\Record Node 101",session=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=1)#
#run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\SOAMMNSOA_2020-12-01_14-09-21_m0003\Record Node 101",session=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=2)
#run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\SOANONSOA_2020-11-16_15-14-52_m0002\Record Node 101",session=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=1)

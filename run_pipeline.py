# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 18:27:13 2021

@author: PC
"""
from mne import set_eeg_reference, events_from_annotations
from mne.epochs import Epochs
from autoreject import AutoReject, Ransac
import numpy as np
from mne.preprocessing.ica import ICA, corrmap, read_ica
import os
from pathlib import Path
from matplotlib import pyplot as plt, patches


def run_pipeline_n1(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea\TESTAUD_2020-11-16_15-14-52_m0002_1sandothers\Record Node 101",experiment=1,out_folder='C:\\Users\PC\Desktop\mea_anal',b=1):
    import mne#this method cycles through all recordings for one recording session. b input is to tell where recordings'number starts from
    global _out_folder
    _out_folder = out_folder
    from glob import glob
    Folder=folder
    experiment=experiment
    Files = sorted(glob(Folder+'/**/*.dat', recursive=True))

    for Fi,File in enumerate(Files):
        b=b#2 for sec animal first series.
        recording=Fi+b#2*second animal
        print('now processing recording'+str(Fi+b))
        indexo,filenamex,soa, stimdur,mean,data,datar,triggers,timestamps,sfreq=loadOep(folder, experiment, recording)#get data
        
        if len(triggers)<100:# (to avoid processing recordings that are botched)
            continue
        raw,events=MNEify(indexo,data,timestamps,triggers,sfreq)
        raw=filtering(raw, notch=50, highpass=None, lowpass=None,
              fir_window="hamming", fir_design="firwin")
        raw=filtering(raw, notch=None, highpass=None, lowpass=100,
              fir_window="hamming", fir_design="firwin")
        reject = dict(
              eeg=280e-6, # V (EEG channels)

              )
        if np.array(indexo).any():
            event_dict = {'stim_standard': 1, 'omission': 0}
            epochs = mne.Epochs(raw, events, event_id=event_dict, tmin= -0.06, reject=reject, baseline=None, preload=True, flat = None, proj=False,reject_by_annotation=False)
            while len(epochs.events)<len(events)*0.8:
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
                epochs = mne.Epochs(raw, events, tmin= -0.06,event_id=event_dict, reject=reject, baseline=None, preload=True, flat = None, proj=False,reject_by_annotation=False, verbose='WARNING')
                data=np.delete(data,todrop[0,0],1)
            epochedstim=np.array(epochs['stim_standard']) 
            epochedomis=np.array(epochs['omission']) 
            
         #   evoked=epochs['omission'].average()
	
         #   myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_evoked-ave.fif'
	
        #    evoked.save(myfile) 
         #   evokede=epochs['omission'].standard_error()
	
         #   myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_evoked-ste_ave.fif'
        #    evokede.save(myfile) 
         #   evoked=epochs['stim_standard'].average()
	#
         #   myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_evoked-ave.fif'
	#
         #   evoked.save(myfile) 
        #    evokede=epochs['stim_standard'].standard_error()
	#
         #   myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_evoked-ste_ave.fif'
         #   evokede.save(myfile) 
            myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_rms_np.npy'

            from pathlib import Path
            sem,std,  bpmean, confidencel,confidenceh=Bootsrapper_snr(epochsdata=epochs['stim_standard']._data)
            stats=[sem,std,  bpmean, confidencel,confidenceh]
            np.save(filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_stats.npy',stats)
            
            sem,std,  bpmean, confidencel,confidenceh=Bootsrapper_snr(epochsdata=epochs['omission']._data)
            stats=[sem,std,  bpmean, confidencel,confidenceh]
            np.save(filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_stats.npy',stats)

        
            my_file = Path(myfile)
            # if my_file.is_file():
            #     sem,std,  bpmean, confidencel,confidenceh=Bootsrapper_snr(epochsdata=epochs['stim_standard']._data)
            #     stats=[sem,std,  bpmean, confidencel,confidenceh]
            #     np.save(filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_stats.npy',stats)
                
            #     sem,std,  bpmean, confidencel,confidenceh=Bootsrapper_snr(epochsdata=epochs['omission']._data)
            #     stats=[sem,std,  bpmean, confidencel,confidenceh]
            #     np.save(filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_stats.npy',stats)


                
            # else:
                
            #     print('save result files and plot stuff')
            #     rmsnp,snr_before_event, snr_later, snr_avg_inverted,peaksnr,peaksnrwide,peakamp,posindex,negindex=snr_peak_by_trial(epochs,epochedstim)#takes lon, save array. also look weird
            #     my_file = Path(myfile)
            #     np.save(myfile, peakamp)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_rms_np.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, rmsnp)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_posindex.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, posindex)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_negindex.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, negindex)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_snr_before_event.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_before_event)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_snr_later.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_later)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_snr_avg_inverted.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_avg_inverted)    
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_peaksnr.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, peaksnr)    
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stan_peaksnrwide.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, peaksnrwide)    
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'indexo.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, np.array(indexo))    

            #     rmsnp,snr_before_event, snr_later, snr_avg_inverted,peaksnr,peaksnrwide,peakamp,posindex,negindex=snr_peak_by_trial(epochs,epochedomis)#takes lon, save array. also look weird
            #     my_file = Path(myfile)
            #     np.save(myfile, peakamp)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_rms_np.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, rmsnp)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_posindex.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, posindex)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_negindex.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, negindex)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_snr_before_event.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_before_event)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_snr_later.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_later)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_snr_avg_inverted.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_avg_inverted)    
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_peaksnr.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, peaksnr)    
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'omis_peaksnrwide.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, peaksnrwide)    



        else:
           
            epochs = mne.Epochs(raw, events, tmin= -0.06, reject=reject, baseline=None, preload=True, flat = None, proj=False,reject_by_annotation=False)
            #dropping=[]
            while len(epochs.events)<len(events)*0.7:
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
	
  #          myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_evoked-ave.fif'
	
    #        evoked.save(myfile) 
       #     evokede=epochs.standard_error()
	#
       #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_evoked-ste_ave.fif'
        #    evokede.save(myfile)     
            myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_rms_np.npy'

            from pathlib import Path
            print('load result files')
            sem,std,  bpmean, confidencel,confidenceh=Bootsrapper_snr(epochsdata=epochs['1']._data)
            stats=[sem,std,  bpmean, confidencel,confidenceh]
            np.save(filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stats.npy',stats)
        
            my_file = Path(myfile)
            # if my_file.is_file():
            #     print('load result files')
            #     sem,std,  bpmean, confidencel,confidenceh=Bootsrapper_snr(epochsdata=epochs['1']._data)
            #     stats=[sem,std,  bpmean, confidencel,confidenceh]
            #     np.save(filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'stats.npy',stats)

            # else:
                
            #     print('save result files and plot stuff')
            #     rmsnp,snr_before_event, snr_later, snr_avg_inverted,peaksnr,peaksnrwide,peakamp,posindex,negindex=snr_peak_by_trial(epochs,epoched)#takes lon, save array. also look weird
            #     my_file = Path(myfile)
            #     np.save(myfile, peakamp)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_rms_np.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, rmsnp)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_posindex.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, posindex)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_negindex.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, negindex)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_snr_before_event.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_before_event)            
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_snr_later.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_later)
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_snr_avg_inverted.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, snr_avg_inverted)    
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_peaksnr.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, peaksnr)    
            #     myfile=filenamex[:-3]+str(round(soa,3))+'_'+str(round(stimdur,3))+'_peaksnrwide.npy'
            #     my_file = Path(myfile)
            #     np.save(myfile, peaksnrwide)    
            
            

def snr_plotter(folder="Z:\\Alessandro_Braga\\n1data from MEA and TDT, first round\mea", outputfolder="C:\\Users\PC\Desktop\snr_plots",which=0, newfreq=500,oldfreq=30000):#which to 0 to select best omis channel
    import matplotlib.pyplot as plt
    from glob import glob
    import mne
    Files = sorted(glob(folder+'/*/*/*/*/*/*\\*.npy', recursive=True))
    u=[]
    for f, file in enumerate(Files):
        if not 'snr' in file:
            u.append(file)
    for name in u:
        Files.remove(name)#now we have tthe list of all the snr files. 
    b=[]
    for f in range(0,len(Files),5):
        print(f)
        if 'stan' in Files[f]:
            continue
        peaksnr=np.load(Files[f])
        peaksnrwide=np.load(Files[f+1])
        prebase_snr=np.load(Files[f+3])
        latebase_snr=np.load(Files[f+4])

        if 'omis' in Files[f]:
            prebaseo_snr=np.load(Files[f+3])
            prebases_snr=np.load(Files[f+3+5])
            condition='omission'
            filevokavg=sorted(glob(Files[f][:Files[f].find('continuous.')]+'*.fif',recursive=True))[0]
            filevokste=sorted(glob(Files[f][:Files[f].find('continuous.')]+'*.fif',recursive=True))[1]
            evokavgo=mne.Evoked(filevokavg, condition=condition, proj=True, kind='average', allow_maxshield=False, verbose=None)
            evoksteo=mne.Evoked(filevokste, condition=condition, proj=True, kind='standard_error', allow_maxshield=False, verbose=None)
            condition='stim_standard'
            filevokavg=sorted(glob(Files[f][:Files[f].find('continuous.')]+'*.fif',recursive=True))[2]
            filevokste=sorted(glob(Files[f][:Files[f].find('continuous.')]+'*.fif',recursive=True))[3]
            evokavgs=mne.Evoked(filevokavg, condition=condition, proj=True, kind='average', allow_maxshield=False, verbose=None)
            evokstes=mne.Evoked(filevokste, condition=condition, proj=True, kind='standard_error', allow_maxshield=False, verbose=None)
            higherpreindo=np.argwhere(prebaseo_snr[:,-1]==np.max(prebaseo_snr[:,-1]))
            lowerpreindo=np.argwhere(prebaseo_snr[:,-1]==np.min(prebaseo_snr[:,-1]))
            higherpreinds=np.argwhere(prebases_snr[:,-1]==np.max(prebases_snr[:,-1]))
            lowerpreinds=np.argwhere(prebases_snr[:,-1]==np.min(prebases_snr[:,-1]))
            avgpreo=prebaseo_snr.mean(0)[20:]        
            avgpres=prebases_snr.mean(0)[20:]   
            ndp=len(evokavgs.data[1])
            if which==0:
                higherpreo=prebaseo_snr[higherpreindo[0]][0][20:]
                higherpres=prebases_snr[higherpreindo[0]][0][20:]
            else:
                higherpreo=prebaseo_snr[higherpreinds[0]][0][20:]
                higherpres=prebases_snr[higherpreinds[0]][0][20:]
        else:
            filevokavg=sorted(glob(Files[f][:Files[f].find('continuous.')]+'*.fif',recursive=True))[0]
            filevokste=sorted(glob(Files[f][:Files[f].find('continuous.')]+'*.fif',recursive=True))[1]
            peaksnr=np.load(Files[f])
            peaksnrwide=np.load(Files[f+1])
            prebase_snr=np.load(Files[f+3])
            latebase_snr=np.load(Files[f+4])
            condition=None
            evokavg=mne.Evoked(filevokavg, condition=condition, proj=True, kind=None, allow_maxshield=False, verbose=None)
            evokste=mne.Evoked(filevokste, condition=condition, proj=True, kind=None, allow_maxshield=False, verbose=None)
            higherpreind=np.argwhere(prebase_snr[:,-1]==np.max(prebase_snr[:,-1]))
            avgpre=prebase_snr.mean(0)[20:]        
            ndp=len(evokavg.data[1])
        
        
            lowerpreind=np.argwhere(prebase_snr[:,-1]==np.min(prebase_snr[:,-1]))


        #problem: some of these channels make no sense. the peak particularly. we will use only rms
        
        #make 1 polt for prerms and one for laterms. start from 20 dp, plot: higher snr, average snr, lower snr with channel names. fit curve.
        #need evoked fif for channel names
        

        
            higherpre=prebase_snr[higherpreind[0]][0][20:]
            lowerpre=prebase_snr[lowerpreind[0]][0][20:]
        color=['#F72585','#B5179E','#7209B7']
        colorfade=['#e382b7','#ad56e0']
        plt.rcParams["figure.figsize"] = [16,8]
        
        rec=Files[f][Files[f].find('recording')+9:Files[f].find('recording')+11]
        rec=rec.replace('\\','')
        import matplotlib
        matplotlib.rcParams.update({'font.size': 15})
        matplotlib.rcParams.keys()
        matplotlib.rcParams.update({'axes.labelweight': 'bold'})
        time= np. arange(0,ndp/newfreq,1/newfreq)
        time=time-time[30]
        num=Files[f][Files[f].find('m0'):Files[f].find('m0')+5]
        uu=sorted(glob(Files[f][:Files[f].find('continuous')]+'*.h5',recursive=True))


        if 'omis' in Files[f]:
            par=Files[f][Files[f].find('continuous.')+11:Files[f].find('omis')]

            if which==0:
                higherpreinds=higherpreindo
                whicht= '(best O)'
            else:
                higherpreindo=higherpreinds
                whicht='(best S)'
            import tables
            ff=tables.open_file(uu[0])
            ntrial=ff.root._v_attrs._0_Setting_reps
            probdev=str(ff.root._v_attrs._0_Setting_probdev)
            plt.figure()
            plt.suptitle('oMMN '+num+' rec'+rec+'. SOA_stim.duration(sec): '+ par+ '. prob deviant='+probdev+', dev trials=' +str(len(prebaseo_snr[1]))+ ', stan trials=' +str(len(prebases_snr[1]))+whicht)
            ax1 = plt.subplot(212)
            ax1.set_title('ERP channel '+evokavgo.ch_names[higherpreindo[0][0]],fontweight="bold")
            evoksteo.data=evoksteo.data*1000000
            evokavgo.data=evokavgo.data*1000000
            evokstes.data=evokstes.data*1000000
            evokavgs.data=evokavgs.data*1000000
            ax1.plot(time, evokavgo.data[higherpreindo][0,0,:],label='omit',linewidth=3, color=color[0])
            plt.fill_between( time, evoksteo.data[higherpreindo][0,0,:]+  evokavgo.data[higherpreindo][0,0,:],   evokavgo.data[higherpreindo][0,0,:]-evoksteo.data[higherpreindo][0,0,:],color=color[0],alpha=0.3,label='STEo')
            ax1.plot(time, evokavgs.data[higherpreinds][0,0,:],label='stand',linewidth=3, color=color[2])
            plt.fill_between( time, evokstes.data[higherpreinds][0,0,:]+  evokavgs.data[higherpreinds][0,0,:],   evokavgs.data[higherpreinds][0,0,:]-evokstes.data[higherpreinds][0,0,:],color=color[2],alpha=0.3,label='STEs')
            plt.xlabel('S',fontweight="bold")
            plt.ylabel('uV',fontweight="bold")
            plt.ylim((-20,+20))
            plt.xlim((-0.06,+0.25))
            plt.legend(loc=4)

            ax2 = plt.subplot(221)
            
            ax2.plot( np.arange(20,20+len(higherpreo),1),higherpreo,color=color[0],linewidth=3,label='ch '+evokavgo.ch_names[higherpreindo[0][0]])
            ax2.plot( np.arange(20,20+len(avgpreo),1),avgpreo,color=colorfade[0],linewidth=3,label='avg all')
            ax2.set_title('SNR omission',fontweight="bold")
            plt.legend()
            plt.ylabel('RMSres/RMSbase',fontweight="bold")
            plt.xlabel('trials',fontweight="bold")   
            ax3 = plt.subplot(222)
            
            ax3.plot( np.arange(20,20+len(higherpres),1),higherpres,color=color[2],linewidth=3,label='ch '+evokavgs.ch_names[higherpreindo[0][0]])
            ax3.plot( np.arange(20,20+len(avgpres),1),avgpres,color=colorfade[1],linewidth=3,label='avg all')
            ax3.set_title('SNR standard',fontweight="bold")
            plt.ylabel('RMSres/RMSbase',fontweight="bold")
            plt.xlabel('trials',fontweight="bold")

            plt.subplots_adjust(hspace = 0.3)



            plt.legend()
            plt.savefig(outputfolder+'//oMMN '+num+whicht+' rec'+rec+'_'+ par+'.png')
            plt.close()
        else:
            par=Files[f][Files[f].find('continuous.')+11:Files[f].find('peak')-1]

            if which==0:
                plt.figure()
                plt.suptitle('N1 '+num+' rec'+rec+'. SOA_stim.duration(sec): '+ par+ '. n trials=' +str(len(prebase_snr[1])))

                ax1 = plt.subplot(211)
                ax1.set_title('ERP channel '+evokavg.ch_names[higherpreind[0][0]],fontweight="bold")
                evokste.data=evokste.data*1000000
                evokavg.data=evokavg.data*1000000

                ax1.plot(time, evokavg.data[higherpreind][0,0,:],label='best ch',linewidth=3, color=color[2])
                plt.fill_between( time, evokste.data[higherpreind][0,0,:]+  evokavg.data[higherpreind][0,0,:],   evokavg.data[higherpreind][0,0,:]-evokste.data[higherpreind][0,0,:],color=color[2],alpha=0.3,label='STEbc')
                
                
                ax1.plot(time, evokavg.data.mean(0),label='avg all',linewidth=3, color=color[1])
                plt.fill_between( time, evokste.data.mean(0)+  evokavg.data.mean(0),   evokavg.data.mean(0)-evokste.data.mean(0),color=color[1],alpha=0.3,label='STEavg')
                plt.xlabel('S',fontweight="bold")
                plt.ylabel('uV',fontweight="bold")
                plt.ylim((-20,+20))
                plt.xlim((-0.06,+0.25))
                plt.legend(loc=4)
                
                ax2 = plt.subplot(212)
                
                ax2.plot( np.arange(20,20+len(higherpre),1),higherpre,color=color[2],linewidth=3,label='ch '+evokavg.ch_names[higherpreind[0][0]])
                ax2.plot( np.arange(20,20+len(avgpre),1),avgpre,color=color[1],linewidth=3,label='avg all')
                ax2.set_title('SNR',fontweight="bold")
                plt.legend()
                plt.ylabel('RMSres/RMSbase',fontweight="bold")
                plt.xlabel('trials',fontweight="bold")
                plt.subplots_adjust(hspace = 0.3)
                    
              
                
    
                plt.show()
                plt.savefig(outputfolder+'//N1 '+num+' rec'+rec+'_'+ par+'.png')
                plt.close()
import sys
import os
import neuron
# sys.path.append('/home/nbrake/pkgs/LFPy-2.2.4')
import LFPy
import numpy as np
from scipy.io import savemat
from importlib import reload

def read_templatename(templatefile):
    """ Load spike times from the file specified
    by the variable spikingFile """
    with open(templatefile, newline='') as file:
         lines = file.readlines()
         for line in lines:
            if(line.find('begintemplate')>-1):
                temp = line.split(' ')
                templatename = temp[1][:-1]
                print(templatename)
                return templatename

    templatename = -1
    return templatename

def main(cellID,EI=5):
    path0 = os.path.join(r'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models',cellID)
    os.chdir(path0)
    os.makedirs(os.path.join(path0,'matlab_recordings'), exist_ok=True)

    path1 = os.path.join(path0,'morphology.hoc')
    templatefile = os.path.join(path0,'template.hoc')
    templatename = read_templatename(templatefile)

    cellParameters = {
        'morphology' : path1,
        'templatefile' :  templatefile,
        'templatename' :  templatename,
        'templateargs' :  0,
        'Ra': 100,
        'celsius': 34,
        'v_init': -65
    }
    eSynParams = {
        "idx": 0,
        "e": 0,
        "syntype": "Exp2Syn",
        "tau1": 0.3,
        "tau2": 1.8,
        "weight": 0.0007,
        "record_current": False
    }
    iSynParams = {
        "idx": 0,
        "e":-80,
        "syntype": "Exp2Syn",
        "tau1": 1,
        "tau2": 10,
        "weight": 0.0007,
        "record_current": False
    }

    neuron.h.load_file("biophysics.hoc")
    cell = LFPy.TemplateCell(**cellParameters)
    cell.tstop = 2000

    L = 0
    for sec in neuron.h.SectionList[0]:
        # Count dendrite length
        idx = sec.name().find('dend')
        if(idx>-1):
            L += sec.L

    # Add excitatory synapses
    mE = int(np.floor(L*1))
    lamE = 1.75/(1+EI*0.15)
    idcsE = cell.get_rand_idx_area_and_distribution_norm(section='allsec',nidx=mE)
    synTimesE = list()
    syn = list()
    for idx in idcsE:
        N = np.random.poisson(lamE*cell.tstop/1000)
        eSynParams['idx'] = idx
        syn.append(LFPy.Synapse(cell, **eSynParams))
        synTimes = cell.tstop*np.random.random(N)
        syn[-1].set_spike_times(synTimes)
        synTimesE.append(synTimes)

    # Add inhibitory synapses
    mI = int(np.floor(L*0.15))
    lamI = EI*lamE
    idcsI = cell.get_rand_idx_area_and_distribution_norm(section='allsec',nidx=mI)
    synTimesI = list()
    for idx in idcsI:
        N = np.random.poisson(lamI*cell.tstop/1000)
        iSynParams['idx'] = idx
        syn.append(LFPy.Synapse(cell, **iSynParams))
        synTimes = cell.tstop*np.random.random(N)
        syn[-1].set_spike_times(synTimes)
        synTimesI.append(synTimes)

    cdm = LFPy.CurrentDipoleMoment(cell=cell)
    cell.simulate(rec_somav=True,probes=[cdm])
    t = cell.tvec.reshape([-1,1])
    v = cell.somav.reshape([-1,1])

    saveFile = os.path.join(path0,'matlab_recordings','synaptic_input_EI' + str(EI).zfill(2) + '.mat')
    savemat(saveFile, {'time':t,'voltage':v,'dipoles':cdm.data.T})

    ###############################
    ########## Resimulate #########
    cell.strip_hoc_objects()
    reload(neuron)

    neuron.h.load_file("biophysics.hoc")
    cell = LFPy.TemplateCell(**cellParameters)
    cell.tstop = 2000
    ###############################

    # Add excitatory synapses
    syn = list()
    for i,idx in enumerate(idcsE):
        eSynParams['idx'] = idx
        syn.append(LFPy.Synapse(cell, **eSynParams))
        syn[-1].set_spike_times(synTimesE[i])

    # Add inhibitory synapses
    for i,idx in enumerate(idcsI):
        iSynParams['idx'] = idx
        syn.append(LFPy.Synapse(cell, **iSynParams))
        syn[-1].set_spike_times(synTimesI[i])

    for sec in neuron.h.SectionList[0]:
        idx = sec.name().find('dend')
        if(idx==-1):
            # Get biophysical mechanisms
            mechs = list()
            mname = neuron.h.ref("")
            mt = neuron.h.MechanismType(0)
            for i in range(int(mt.count())):
                mt.select(i)
                mt.selected(mname)
                name = mname[0]
                if name in dir(sec(0.5)):
                    if not name.endswith("_ion"):
                        mechs.append(name)
            # set Na conductaqnces to 0
            Na_idcs = np.argwhere([x.find('Na')>=0 for x in mechs]).flatten()
            for idx in Na_idcs:
                Na_mech = mechs[idx]
                par_name = 'g'+Na_mech+'bar_'+Na_mech
                setattr(sec,par_name,0)

    cdm = LFPy.CurrentDipoleMoment(cell=cell)
    cell.simulate(rec_somav=True,probes=[cdm])
    t = cell.tvec.reshape([-1,1])
    v = cell.somav.reshape([-1,1])

    saveFile = os.path.join(path0,'matlab_recordings','synaptic_input_EI' + str(EI).zfill(2) + '_passive.mat')
    savemat(saveFile, {'time':t,'voltage':v,'dipoles':cdm.data.T})

if __name__ == '__main__':
    if(len(sys.argv)==2):
        path0 = os.path.join(r'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models',sys.argv[1])
        file = os.path.join(path0,'EI_ratio.csv')
        with open(file,"r") as f:
            EI = eval(f.readline())
    elif(len(sys.argv)==3):
        EI_vec = [1,1.5,2.1,3.1,4.5,6.6,9.7,14.1,20.6,30]
        idx = int(sys.argv[2])-1
        EI = EI_vec[idx]

    main(sys.argv[1],EI)
import sys
import os
import neuron
import LFPy
import numpy as np

cellID = 'L5_TTPC2_cADpyr232_1'

path0 = os.path.join(r'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\other_code\neuron_models',cellID)
mech_folder = os.path.join(path0,'mechanisms')
os.chdir(path0)
os.system("nrnivmodl " + mech_folder)
os.system("mknrndll " + mech_folder)
os.makedirs(os.path.join(path0,'matlab_recording'), exist_ok=True)

path1 = os.path.join(path0,'morphology.hoc')
templatefile = os.path.join(path0,'template.hoc')

cellParameters = {
    'morphology' : path1,
    'templatefile' :  templatefile,
    'templatename' :  'cADpyr232_L5_TTPC2_8052133265',
    'templateargs' :  0,
    'Ra': 100,
    'passive' : False,
    'celsius': 34,
    'v_init': -65
}
neuron.h.load_file("biophysics.hoc")

cell = LFPy.TemplateCell(**cellParameters)
cell.tstop = 3000

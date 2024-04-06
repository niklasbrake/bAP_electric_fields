import sys
sys.path.append('/lustre04/scratch/nbrake/code/simulation_code')
import os
import neuron
sys.path.append('/home/nbrake/pkgs/LFPy-2.2.4')
sys.path.append('/home/nbrake/pkgs/hoc2swc-master')
import LFPy
import numpy as np
from scipy.io import savemat
from hoc2swc import neuron2swc

# nrnivmodl mechanisms

if __name__ == '__main__':
    # cellID = 'L23_DBC_bIR215_4'
    cellID = sys.argv[1]
    print(cellID)

    path0 = os.path.join('/lustre04/scratch/nbrake/data/simulations/unitary_AP',cellID)
    os.chdir(path0)
    os.makedirs(os.path.join(path0,'matlab_recordings'), exist_ok=True)

    path1 = os.path.join(path0,'morphology.hoc')
    templatefile = os.path.join(path0,'template.hoc')


    with open(templatefile, newline='') as file:
         lines = file.readlines()
         for line in lines:
            if(line.find('begintemplate')>-1):
                temp = line.split(' ')
                templatename = temp[1][:-1]
                print(templatename)

    cellParameters = {
        'morphology' : path1,
        'templatefile' :  templatefile,
        'templatename' :  templatename,
        'templateargs' :  0,
        'Ra': 100,
        'passive' : False,
        'celsius': 34,
        'v_init': -65
    }
    neuron.h.load_file("biophysics.hoc")

    cell = LFPy.TemplateCell(**cellParameters)
    cell.tstop = 3000

    neuron2swc("out.swc")

    import getMorphoSegments
    pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs("out.swc")
    S = np.array(segs,dtype=object)
    savemat('morphData.mat', {'connections':connections,'data':morphData,'segs':S})


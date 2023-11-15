from pathlib import Path
from pyLCIO import IOIMPL, IMPL, UTIL, EVENT

import os
import numpy as np
from pattern.detector_hit import DetectorHit
from pattern.multiplet import Multiplet
 

file = '/nfs/dust/luxe/user/spatarod/key4hepcsv/p0xi3BX0000/' \
       'event000000000-simhits-LUXE_example_geant4/_eigensolver_7q_impact_list_reverse'
    
reco_raw = 'reco_xplet_list.npy'
reco_ambiguity_solved = 'reco_xplet_list_ambiguity_solved.npy'

wrt = IOIMPL.LCFactory.getInstance().createLCWriter()
wrt.open(f"{file}/ACTS", EVENT.LCIO.WRITE_NEW)

string_encoding = 'system:1,side:1,layer:3,module:5,sensor:0,x:32:-16,y:-16'
xplets = np.load(f'{file}/{reco_raw}',allow_pickle=True)

for i, xplet in enumerate(xplets):
    col = IMPL.LCCollectionVec(EVENT.LCIO.TRACKERHITPLANE)
    evt = IMPL.LCEventImpl()
    for x,y,z,cell_id in zip(xplet.x, xplet.y, xplet.z, xplet.cell_id):
        hit = IMPL.TrackerHitPlaneImpl()
        hit.setdU(0.005)
        hit.setdV(0.005)
        hit.setPosition(np.array([1e3 * x, 1e3 * y, 1e3 * z]))
        hit.setCellID0(cell_id)
        param = col.parameters()
        param.setValue(EVENT.LCIO.CellIDEncoding, string_encoding)
        col.addElement(hit)
            
    evt.setEventNumber(i)
    evt.addCollection(col, 'SiTrackerHits')
    wrt.writeEvent(evt)
    
wrt.close()
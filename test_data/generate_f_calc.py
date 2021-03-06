from __future__ import division
import mmtbx.command_line.fmodel
import mmtbx.utils
from iotbx import file_reader
import math
from cctbx import miller
from libtbx import easy_pickle
import os
import random
import sys

args = sys.argv[1:]
n_frames=int(args[0])
resulotion=float(args[1])

pdb_in = file_reader.any_file("1JB0.pdb", force_type="pdb")
pdb_in.assert_file_type("pdb")
xray_structure = pdb_in.file_object.xray_structure_simple()
xray_structure.show_summary()
phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
params2 = phil2.extract()
# adjust the cutoff of the generated intensities to assure that
# statistics will be reported to the desired high-resolution limit
# even if the observed unit cell differs slightly from the reference.
params2.high_resolution = resulotion
params2.output.type = "real"
f_model = mmtbx.utils.fmodel_from_xray_structure(
  xray_structure = xray_structure,
  f_obs          = None,
  add_sigmas     = True,
  params         = params2).f_model
i_model = f_model.as_intensity_array().map_to_asu()
#print len(i_model.indices())
count=0
write_I=open("I_reference_"+str(n_frames)+".db","w+")
for i, j in zip(i_model.indices(),i_model.data()):
 print>>write_I, j
 count=count+1
print "Total number of miller indecis: ", count
write_I.close()
random.seed(0)
frame=0
write_G=open("G_reference_"+str(n_frames)+".db","w+")
write=open("PSI_simulated_observations_"+str(n_frames)+".db","w+")
for files in os.listdir("/reg/neh/home/mamin03/PSI/integration"):
 print "Frame number ", frame
 if frame==n_frames: break
 d=easy_pickle.load("/reg/neh/home/mamin03/PSI/integration/"+files)
 A=d['observations'][0].as_non_anomalous_array().map_to_asu()
 matches = miller.match_multi_indices(
           miller_indices_unique=i_model.indices(),
           miller_indices=A.indices())
 scale=random.random()*10 #scale factor for each frame (simulated G)
 print>>write_G, scale
 for pair in matches.pairs():
#  print pair[0] , i_model.indices()[pair[0]], A.indices()[pair[1]], scale*i_model.data()[pair[0]], i_model.sigmas()[pair[0]]
  scaled_I=scale*i_model.data()[pair[0]]
  sigma=math.sqrt(scaled_I)
  #scaled_I=scaled_I+random.normalvariate(0, sigma)
  print>> write, pair[0], scaled_I, sigma, 0, 0, frame
 frame=frame+1
write.close()
write_G.close()

datafile="PSI_simulated_observations_"+str(n_frames)+".db"
with open(datafile) as f:
    file_sorted = sorted((x for x in f),
                         key=lambda z:(int(z.split(" ")[0]),int(z.split(" ")[5])),
                         reverse=False)
f.close()
w=open(datafile, 'w+')
for i in file_sorted:w.write(i)
w.close()

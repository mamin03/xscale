from __future__ import division
from cxi_xdr_xes.XScale import Load_object, read_objects
from read_files import read_cxi_merge, correlation
import sys
from cctbx.array_family import flex
import time
import iotbx.phil

def get_Ih(hkl_index, I, W, frames, G, n_Miller):
  index=hkl_index[0]
  Nsum=0
  Dsum=0
  Ih=flex.double(n_Miller)
  for i in range(len(I)):
    if index==hkl_index[i]:
      Nsum=Nsum+I[i]*W[i]*G[frames[i]]
      Dsum=Dsum+W[i]*G[frames[i]]*G[frames[i]]
    else:
      Ih[index]=Nsum/Dsum
      index=hkl_index[i]
      Nsum=I[i]*W[i]*G[frames[i]]
      Dsum=W[i]*G[frames[i]]*G[frames[i]]
  Ih[-1]=Nsum/Dsum
  return Ih


master_phil_str = """
xscale {
  method = *G_only G_cycles G_wI G_wI_LM
    .type = choice
    .help = "The method for scaling."
  eps = 1.e2
    .type = float
    .help = convergence threshold.
  n_cycles = 10
    .type = int
    .help = number of minimization cycles (only for the G_cycles method).
  ref_G = None
    .type = str
    .help = the file contains real scale factors
 ref_I = None
    .type = str
    .help = the file contains real intensities
  plot = *No Yes
    .type = choice
    .help = ploting option
}
"""

master_phil = iotbx.phil.parse(master_phil_str)
def run(args):
  processed = iotbx.phil.process_command_line(
    args=args, master_string=master_phil_str)
  args = processed.remaining_args
  work_params = processed.work.extract().xscale
  processed.work.show()
  assert len(args) == 1
  observations=args[0]
  choose_method=work_params.method
  assert choose_method in ("G_only", "G_cycles", "G_wI", "G_wI_LM")
  if choose_method=="G_only":   method=1
  if choose_method=="G_cycles": method=2
  if choose_method=="G_wI":     method=3
  if choose_method=="G_wI_LM":  method=4

  eps    = work_params.eps
  ref_G  = work_params.ref_G
  ref_I  = work_params.ref_I
  cycles = work_params.n_cycles
  mi=None
  f=open(observations)
  (I, W, H, Miller_list, F)=read_cxi_merge(f)
  t1=time.time()
  if method==1:
    from cxi_xdr_xes.XScale.minimizers.scaling_minimizer_G import xscale
    print "Loading objects......."
    n_frames = max(F)+1
    (MI, Fr) = Load_object(I, W, H, Miller_list, F, n_frames)
    fit=xscale(MI, Fr, n_frames, eps)
    Gm=fit.x
    (Ih,mi)=read_objects(MI)

  if method==2:
    from cxi_xdr_xes.XScale.minimizers.scaling_minimizer_G_cycles import xscale
    N_miller=max(H)+1
    Gm=flex.double(max(F)+1,1)
    t1=time.time()
    Ih=get_Ih(H, I, W, F, Gm, N_miller)
    for i in range(cycles):
      fit=xscale(I, W, H, F, flex.double(Gm), flex.double(Ih), eps)
      Gm=fit.x
      Ih=get_Ih(H, I, W, F, Gm, N_miller)

  if method==3:
    from cxi_xdr_xes.XScale.minimizers.scaling_minimizer_wI_G import xscale
    N_miller=max(H)+1
    G=flex.double(max(F)+1,1)
    Ih=get_Ih(H, I, W, F, G, N_miller)
    fit=xscale(I, W, H, F, flex.double(G), flex.double(Ih), eps)
    Ih=fit.x[len(G):]*flex.double(Ih)
    Gm=fit.x[:len(G)]

  if method==4:
    from scitbx.lstbx import normal_eqns_solving
    from cxi_xdr_xes.XScale.minimizers.scaling_LM_minimizer_wI_G import xscale
    N_miller=max(H)+1
    G=flex.double(max(F)+1,1)
    Ih=get_Ih(H, I, W, F, G, N_miller)
    fit=xscale(I, W, H, F, G, Ih)
    iterations = normal_eqns_solving.levenberg_marquardt_iterations(
          fit,
          gradient_threshold = eps)
    Ih=fit.x[len(G):]*flex.double(Ih)
    Gm=fit.x[:len(G)]
  print "Total time", time.time()-t1
  if ref_G or ref_I:
    correlation(ref_G, ref_I, Gm, Ih, mi, False)
if __name__ == '__main__':
  run(sys.argv[1:])

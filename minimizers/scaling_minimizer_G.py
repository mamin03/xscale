from __future__ import division
from cctbx.array_family import flex
from cxi_xdr_xes.XScale import get_f_g
from scitbx import lbfgs
class xscale:
  def __init__(self, MI, Fr, n_frames, eps):
    self.counter=0
    self.MI = MI
    self.Fr = Fr
    self.n_frames=n_frames
    self.x = flex.double(self.n_frames,1)
    termination_params = lbfgs.termination_parameters(
                      traditional_convergence_test=True,
                       traditional_convergence_test_eps=eps)
    self.minimizer = lbfgs.run(target_evaluator=self, termination_params=termination_params)
    print "End of minimization: Converged"
  def compute_functional_and_gradients(self):
   (f, g, c) = get_f_g(self.MI, self.x, self.Fr) #C++ code that returns the
   self.c = c
   print 'The value of F', f
   self.counter = self.counter+1
   return f, g


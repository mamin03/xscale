from __future__ import division
from scitbx.array_family import flex
from scitbx import lbfgs

class xscale:
  def __init__(self, I, weight, hkl, frames, G, Ih, eps):
    self.counter=0
    self.I = flex.double(I)
    self.weight = flex.double(weight)
    self.frames = flex.size_t(frames)
    self.hkl = flex.size_t(hkl)
    self.Ih = flex.double(Ih)
    self.x = flex.double(G)
    termination_params = lbfgs.termination_parameters(
                      traditional_convergence_test=True,
                       traditional_convergence_test_eps=eps)
    self.minimizer = lbfgs.run(target_evaluator=self, termination_params=termination_params)
    print "End of minimization: Converged after", self.counter, "steps"

  def compute_functional_and_gradients(self):

   Ih = self.Ih.select(self.hkl)
   Gm = self.x.select(self.frames)
   t = self.I - Gm * Ih
   f = self.weight * t * t
   gterms = -2. * self.weight * self.Ih.select(self.hkl) * t
   g = flex.double(len(self.x),0.)
   g.add_selected(self.frames,gterms)
   f=flex.sum(f)
   print 'value of F is (refining G only):', f, "Trial Number: ", self.counter
   self.counter = self.counter+1
   return f, g

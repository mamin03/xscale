from __future__ import division
from scitbx.array_family import flex
from scitbx import lbfgs

class xscale:
  def __init__(self, I, weight, hkl, frames, G, Ih, eps):
    self.counter=0
    self.n_frames=max(frames)+1
    self.I = flex.double(I)
    self.weight = flex.double(weight)
    self.frames = flex.size_t(frames)
    self.hkl = flex.size_t(hkl)
    self.s = flex.double(Ih)
    self.x = flex.double((len(Ih)+len(G)),1)
    self.function=0
    termination_params = lbfgs.termination_parameters(
                      traditional_convergence_test=True, traditional_convergence_test_eps=eps)
    self.minimizer = lbfgs.run(target_evaluator=self, termination_params=termination_params)
    print "End of minimization: Converged after", self.counter, "steps"

  def compute_functional_and_gradients(self):
    Gm=self.x[:self.n_frames]
    Ih=self.x[self.n_frames:]
    Ih = Ih.select(self.hkl)
    Gm = Gm.select(self.frames)
    Sh = self.s.select(self.hkl)
    t = self.I - (Gm * Ih * Sh)
    f = self.weight * t * t
    gterms_gr = -2. * self.weight * Ih * Sh * t
    Ihterms_gr = -2. * self.weight * Gm * Sh * t
    g = flex.double(len(self.x),0.)
    gGm = flex.double(self.n_frames,0.)
    gIh = flex.double(len(self.s),0.)
    gGm.add_selected(self.frames,gterms_gr)
    gIh.add_selected(self.hkl, Ihterms_gr)
    g=gGm.concatenate(gIh)
    self.counter = self.counter+1
    self.function=flex.sum(f)
    print 'value of F is:', self.function, "Trial Number: ", self.counter
    return self.function, g


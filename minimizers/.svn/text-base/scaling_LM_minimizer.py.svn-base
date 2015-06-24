from __future__ import division
from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns

class xscale(
  normal_eqns.non_linear_ls,
  normal_eqns.non_linear_ls_mixin):
  def __init__(self, I, sigma, hkl, frames, G, Ih):
    self.I = flex.double(I)
    self.sigma = flex.double(sigma)
    self.hkl = flex.size_t(hkl)
    self.frames = flex.size_t(frames)
    self.G = G
    self.Ih = flex.double(Ih)
    self.n_obs= len(I)
    self.n_prm = len(G)
    self.count=0
    super(xscale, self).__init__(n_parameters=self.n_prm)

    self.restart()

  def restart(self):
    self.x = self.G.deep_copy()
    self.old_x = None
    print "Restart"
  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):
    print "Trial number ...", self.count
    residuals = self.sigma*(self.I-(self.x.select(self.frames)*self.Ih.select(self.hkl)))
    self.reset()
    if objective_only:
      self.add_residuals(residuals, weights=None)
    else:
      jacobian = flex.double(flex.grid(self.n_obs, self.n_prm))
      for row in range(self.n_obs):
        for column in range(self.n_prm):
          if column==self.frames[row]:
            h=self.hkl[row]
            jacobian[(row, column)]=-self.sigma[row]*self.Ih[h]
            break
      self.add_equations(residuals, jacobian, weights=None)
    self.count+=1
  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None

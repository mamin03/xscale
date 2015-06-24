from __future__ import division
from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns
from cxi_xdr_xes import calculate_residual_jacobian_wI_G
from scitbx import sparse
import time
class xscale(
 # normal_eqns.non_linear_ls_with_separable_scale_factor,
  normal_eqns.non_linear_ls,
  normal_eqns.non_linear_ls_mixin):
  def __init__(self, I, w, hkl, frames, G, Ih):
    self.I = flex.double(I)
    self.w = flex.double(w)
    self.hkl = flex.size_t(hkl)
    self.frames = flex.size_t(frames)
    self.G = G
    self.Ih = flex.double(Ih)
    self.S = flex.double(Ih)
    self.n_obs= len(I)
    self.n_frames = len(G)
    self.n_prm = self.n_frames + len(self.S)
    self.count=0
    self.weight = flex.double(self.n_obs,1)
    self.jacobian_sparse = sparse.matrix(self.n_obs, self.n_prm)
    super(xscale, self).__init__(n_parameters=self.n_prm)
    self.restart()

  def restart(self):
    self.x = flex.double(self.n_prm,1)
    self.old_x = None
    print "Restart"
  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):
    print "Trial number ...", self.count
    G = self.x[:self.n_frames]
    Ih = self.x[self.n_frames:]
    print "calculating residuales and Jacobians......"
    residuals=flex.double(self.n_obs)
    for  i in range(self.n_obs):
      m=self.frames[i]
      h = self.hkl[i]
      residuals[i] = self.w[i]*(self.I[i] - (G[m] * Ih[h] * self.S[h]));
      self.jacobian_sparse[i, m] = -self.w[i] * Ih[h] * self.S[h];
      self.jacobian_sparse[i, h+self.n_frames] = -self.w[i] * G[m] * self.S[h];
    print "Done .."
    self.reset()
    if objective_only:
      self.add_residuals(residuals, self.weight)
    else:
      print "add equation"
      t1=time.time()
      self.add_equations(residuals, self.jacobian_sparse, self.weight)
      print "time", time.time()-t1
      print "Done add equation"
    self.count+=1
  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None

from __future__ import division
from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns
from cxi_xdr_xes.XScale import calculate_residual_jacobian_wI_G
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
    self.jacobian = flex.double(flex.grid(self.n_obs, self.n_prm))
    super(xscale, self).__init__(n_parameters=self.n_prm)
    self.restart()
    self.build_up()
  def restart(self):
    self.x = flex.double(self.n_prm,1)
    self.old_x = None
    print "Restart"
  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):
    print "Trial number ...", self.count
    Gm = self.x[:self.n_frames]
    Ih = self.x[self.n_frames:]
    print "calculating residuals and Jacobians......"
    (residuals, jacobian) = calculate_residual_jacobian_wI_G(self.I,
                                                             self.w,
                                                             self.hkl,
                                                             self.frames,
                                                             self.S,
                                                             Gm,
                                                             Ih,
                                                             self.n_frames,
                                                             self.jacobian)

#    J=jacobian.as_numpy_array()
#    JT=jacobian.matrix_transpose().as_numpy_array()
#    A=jacobian.matrix_transpose_multiply(jacobian)
#    nonzero=0
#    all=0
#    B=A.as_numpy_array()
#    print B
#    print JT
#    print J
#    for i in A:
#      if i!=0:
#        nonzero+=1
#      all+=1
#    print "Density", nonzero/all
    print "Done .."
    self.reset()
    if objective_only:
      self.add_residuals(residuals, weights=None)
    else:
      self.add_equations(residuals, jacobian, weights=None)
      print "Done add equation"
    self.count+=1
  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None

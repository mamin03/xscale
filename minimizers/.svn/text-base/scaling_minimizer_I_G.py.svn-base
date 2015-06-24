from __future__ import division
from cctbx.array_family import flex
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in


class XScale(lbfgs_with_curvatures_mix_in):
  def __init__(self, MI):
    self.MI = MI
    self.n_frames=self.get_Frames()
    self.x = flex.double(self.n_frames + len(MI))
    for i in range(self.n_frames):
      self.x[i] = 1 # initiate all scale factors to 1
    for i in range(len(MI)):
      self.x[i+self.n_frames] = MI[i].get_avg()  # initiate all intensities to average
    super(XScale, self).__init__(use_curvatures=True)

  def compute_functional_and_gradients(self):
   f = 0
   g = flex.double(len(self.x))
   self.c = flex.double(len(self.x))
   for i in self.MI: # the first summation in Eq.1 (loop over all miller indices)
     for j in i.I_list: # the second summation in Eq.1 (loop over all obsevations)
       m = j.Frame.ID
       j.Frame.G = self.x[m]
       h = i.index
       i.refined_I = self.x[self.n_frames+h]
       t = j.value - j.Frame.G*i.refined_I
       f = f + t**2
       g[m] = g[m] - 2 * i.refined_I * t
       g[self.n_frames + h] = g[self.n_frames + h]-2.0 * j.Frame.G * t;
       self.c[m] = self.c[m] + 2 * i.refined_I * i.refined_I          # calculate the second derivatives here
       self.c[self.n_frames + h] = self.c[self.n_frames + h] + 2.0 * j.Frame.G * j.Frame.G;
   print f
   return f, g

  def curvatures(self):
   print 'Curvature Called'
   return self.c

  def get_Frames(self):
    self.Frames=[]
    for i in self.MI:
       for j in i.I_list:
         exist=0
         if len(self.Frames)==0:
             self.Frames.append(j.Frame)
         else:
           for k in self.Frames:
              if k.ID==j.Frame.ID:
                exist=1
           if exist==0:
              self.Frames.append(j.Frame)
    return len(self.Frames)

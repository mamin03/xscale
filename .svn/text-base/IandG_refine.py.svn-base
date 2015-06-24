from __future__ import division
def get_nframes(F):
  frames=[]
  for i in F:
    if i in frames:
      pass
    else: frames.append(i)
  return len(frames)

def get_avg_I(I,H):
  miller=[]
  for i in H:
    if i in miller:
      pass
    else: miller.append(i)
  I_avg=[]
  for i in miller: I_avg.append(0)
  for i in range(len(miller)):
    count=0
    for j in range(len(H)):
      if H[j]==miller[i]:
       count = count+1
       I_avg[i] = I_avg[i] + I[j]
    if count!=0:
     I_avg[i]=I_avg[i]/count
  return I_avg

def get_miller(H):
  miller=[]
  for i in H:
    if i in miller:
      pass
    else: miller.append(i)
  return len(miller)

from cctbx.array_family import flex
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in


class XScale(lbfgs_with_curvatures_mix_in):
  def __init__(self, I, sig, hkl, frame_ID):
    # Initialise all per-frame scale factors to one.
    self.I = I
    self.sig = sig
    self.hkl = hkl
    self.frame_ID = frame_ID
    self.n_frames = get_nframes(frame_ID)
    self.n_miller = get_miller(hkl)
    self.x = flex.double(self.n_frames + self.n_miller)
    I_avg=get_avg_I(I,H)
#    I_avg=[1000, 1200, 1300, 900]
    for i in range(self.n_frames):
      self.x[i] = 2
    for i in range(self.n_miller):
        self.x[self.n_frames + i] = I_avg[i]
    super(XScale, self).__init__(use_curvatures=True)

  def compute_functional_and_gradients(self):
   f = 0
   g = flex.double(len(self.x))
   for i in range(len(self.I)):
     m = self.frame_ID[i]
     G = self.x[m]
     h = self.hkl[i]
     I_m = self.x[self.n_frames+h]
     t = self.I[i]-G*I_m
     f = f + t**2
     g[m] = g[m] - 2 * I_m * t
     g[self.n_frames + h] = g[self.n_frames + h]-2.0 * G * t;
   print f
   return f, g

  def curvatures(self):
   c = flex.double(len(self.x))
   for i in range(len(self.I)):
     m = self.frame_ID[i]
     G = self.x[m]
     h = self.hkl[i]
     I_m = self.x[self.n_frames+h]
     c[m] = c[m] + 2 * I_m * I_m
     c[self.n_frames + h] = c[self.n_frames + h] + 2.0 * G * G;
   return c

# refactor to XSCALE and test
# keep all tests
# random numbers seeds (set up to 0)
# check in the code Litex
# not enough comments
# matplot is an option
# correlation coeff. as test critrea ("pass" or "fail")
# standard test of realistic data (lys. diffracting to 2a) 1000 Frames.
# C++
# weighting by sigma or varience
# [G+I] or [G]: choose either one with a phil parameter
# correlation coefficent of Imodel using G+I vs G.
# using Tara's simulator to simulate crystals "realistically" in random orientations

# future path:
# propagation of errors.
# B-factor
# test different partiality models
# mosaic block size + mosaicity
# "postrefinement"--refinement of unit cell & orientation.
# taking wavelength dispersion into account.


f=open('obstest.txt','r')
#f=open('test.txt','r')
I=[]
S=[]
H=[]
F=[]
for i in f:
  if i.split(" ")[0]!='I':
    I.append(float(i.split(" ")[0].strip()))
    S.append(float(i.split(" ")[1].strip()))
    H.append(int(i.split(" ")[2].strip()))
    F.append(int(i.split(" ")[3].strip()))
print '---------------------------------------'
print "Original Data"
print '---------------------------------------'
print 'I  ', 'S', 'H', 'F'
for i in range(len(I)):
   print I[i], S[i], H[i], F[i]
I_avg=get_avg_I(I,H)
print '---------------------------------------'
print 'Average I'
print '---------------------------------------'
for i in I_avg: print i
print '---------------------------------------'
print 'Target Function Values'
print '---------------------------------------'
fit = XScale(I,S,H,F)
print '---------------------------------------'
print "Scale Factor G and Scaled I"
print '---------------------------------------'
n_frames = get_nframes(F)
n_miller = get_miller(H)
Imodel=[]
from matplotlib import pyplot as plt
for i in range(n_frames):
  print 'For Frame', i, 'G is ', fit.x[i]
for i in range(n_miller):
  Imodel.append(fit.x[n_frames+i])
  print 'For index', H[i], 'the I is ', fit.x[n_frames+i]

plt.plot(I_avg,Imodel, 'r.')
plt.show()

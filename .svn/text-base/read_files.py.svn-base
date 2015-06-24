from __future__ import division
#from matplotlib import pyplot as plt
from cctbx.array_family import flex

def correlation(G_ref, I_ref, G_model, I_model, mi=None, plot=False):
  if G_ref:
    ref=open(G_ref)
    I_ref_sub=[]
    for i in ref:
      I_ref_sub.append(float(i))
    R=flex.linear_correlation(flex.double(G_model), flex.double(I_ref_sub)).coefficient()
    print "The Gm correlation: ", round(R, 5)
  if I_ref:
    ref=open(I_ref)
    I_ref_sub=[]
    I_ref_match=[]
    for i in ref:
      I_ref_sub.append(float(i))
  if mi:
    for j in  mi:
      I_ref_match.append(I_ref_sub[j])
  else:
    model_match=[]
    for i,j in zip(I_model, I_ref_sub):
     if i!=0:
       model_match.append(i)
       I_ref_match.append(j)
    I_model=model_match
  R=flex.linear_correlation(flex.double(I_model), flex.double(I_ref_match)).coefficient()
  print "The Ih correlation: ", round(R, 5)

#  if(plot==True):
#    plt.plot(I_model, I_ref_sub, 'r.')
#    plt.show()



def read_cxi_merge(f):
  I=[]
  W=[]
  H=[]
  F=[]
  Miller_list=[]

  print 'Reading Observation File..........'
  for i in f:
     token=i.split(" ")
 #   if (float(token[1]))>0: if you want to ignore negative intensities
     I.append(float(token[1]))
     W.append(1/((float(token[2]))*(float(token[2]))))
     H.append(int(token[0]))
     F.append(int(token[5]))
  Miller_list=[x for x in range(max(H)+1)]
  return (I, W, H, Miller_list, F)

def read_random_orientation():
  print "Reading intensities from random orientation"
  Im=[]
  m=open('Im.db')
  for i in m:
    Im.append(float(i))
  f=open('random_orientations.db','r')
  I=[]
  W=[]
  H=[]
  F=[]
  for i in f:
    if i.split(" ")[0]!='I':
      token=i.split(" ")
      I.append(float(token[0]))
      w.append(float(token[1]))
      H.append(int(token[2]))
      F.append(int(token[3]))
  Miller_list=[x for x in xrange(196)]
  return (I, W, H, Miller_list, F, Im)

from __future__ import division

class miller(object):
  """This Class represents a single miller index. Each miller index should carry a list of intensities (I_list)
     from different frames. The function add_I will add new intensity to the list of intensities.
     The Function get_avg returns the the average of all observations for this miller index. This value is
     a starting value if Ih is being refined"""

  def __init__(self, index):  # Initiated with the miller index
    self.index=index
    self.I_list=[]            # Initial an array that carries the Intensties
    self.refined_I=0

  def add_I(self,I):          # To add an Intensity to the list of Intensites of a given miller index
    self.I_list.append(I)

  def get_avg(self):          # Return an average of the Intensties for this miller index
    sum=0
    for i in self.I_list:
       sum=sum+i.value
    return sum/len(self.I_list)

  def Ih(self):
    numerator = 0
    denominator = 0
    for i in self.I_list:
       numerator   = numerator + i.value * i.sigma * i.Frame.G     # Eqn (3) J.Appl.Cryst. (1998). 31, 708-717
       denominator = denominator + i.sigma * i.Frame.G * i.Frame.G
    Ih =  numerator/denominator
    sigma = 1/denominator
    return Ih

  def dIhdG(self, m):               # Drivatives of Eqn (3) using the formula of the rational functions
    numerator = 0
    denominator = 0
    d_numerator = 0
    d_denominator = 0
    d_2_numerator = 0
    d_2_denominator = 0
    for i in self.I_list:
      numerator   = numerator + i.value * i.sigma * i.Frame.G
      denominator = denominator + i.sigma * i.Frame.G * i.Frame.G
      if i.Frame.ID == m:
        d_numerator = i.value * i.sigma
        d_denominator = 2 * i.sigma * i.Frame.G
        d_2_numerator = 0
        d_2_denominator = 2 * i.sigma
    dIhdG = (denominator * d_numerator - numerator * d_denominator) / (denominator * denominator)
    dIh2dG = ((denominator * denominator)*(d_numerator*d_denominator - (numerator * d_2_denominator+d_denominator*d_numerator))-(2 * denominator * d_denominator * (denominator * d_numerator - numerator * d_denominator))) / (denominator * denominator)**2
    return dIhdG, dIh2dG



class Frame:

  def __init__(self, ID): # Each frame has an ID and Scale Factor initiated to 1
    self.ID = ID
    self.G  = 1
    self.B  = 0  # for the B factor (not in the code yet)

  def set_Scale_Factor(self, Scale_Factor): # you will need to set the scale factor after minimization
    self.G = Scale_Factor


class Intensity:     # we will need to add a partialities and wavelengh (for the 2 color)

  def __init__(self, value, sigma, Frame):  # Initiated by the value of I and sigma and Frame object.
    self.value = value
    self.sigma = sigma
    self.Frame = Frame
    self.theta = 0     # to calculate the lorentz and polarization
    self.lamda = 0     # scaling for 2 colors

  def get_Ls():         # will return the lorentz correction
     return

  def get_Lp():         # will return polarization correction
     return

  def weight_I(self, weight): # you may need to wieght the intensities
    self.value_weighted = self.value * weight

  def scale_I(self, Frame): # you may need to scale the intensities
    self.value_scaled = self.value * Frame.G

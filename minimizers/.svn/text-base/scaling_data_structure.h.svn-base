#include <vector>
#include <iostream>
struct derivatives {

double first_dirv, second_dirv;
derivatives( double x,  double y ) : first_dirv(x), second_dirv(y) {}
};

class Miller;

class Frame
{
 public:
    int ID;
    double G, B;
    std::vector< Miller*> miller_list;
    bool exist(Miller *m) {
      for (int i=0; i<miller_list.size(); i++) {
        if(m==miller_list[i])
        return true;
      }
    return false;
    }
};


class Intensity
{

 public:
    double value;
    double sigma;
    Frame  *frame;
};

class Miller
{
    /*This Class represents a single miller index. Each miller index should carry a list of intensities (I_list)
     from different frames. The function add_I will add new intensity to the list of intensities.
     The Function get_avg returns the the average of all observations for this miller index. This value is
     a starting value if Ih is being refined*/
  public:
    int index;
    double Im;
    double numerator;
    double denominator;
    std::vector< Intensity > I_list;
    void add_I (Intensity I)
    {
     I_list.push_back(I);
    }

    double get_avg() //needed for refining I.
    {               //he initial guess is the average of intensities for a given miller index
     double sum=0;
     for (int i=0; i<I_list.size(); i++)
     {
      sum=sum+I_list[i].value;
     }
    return sum/I_list.size();
    }
    double Ih()  //return the value of Ih based on equation 3
    {
     double SigmaIh=0;
     numerator=0;
     denominator = 0;
     for (int i=0; i<I_list.size(); i++)
     {
      numerator = numerator + I_list[i].value * I_list[i].sigma * I_list[i].frame->G;
      denominator = denominator + I_list[i].sigma * I_list[i].frame->G * I_list[i].frame->G;
     }
     Im = numerator/denominator;
     SigmaIh = 1/denominator;
     return Im;
    }
    derivatives dIhdG(int m) //return the first and the second dirvatives of Ih
    {
    double d_numerator = 0;
    double d_denominator = 0;
    double d_2_numerator = 0;
    double d_2_denominator = 0;
    double dIhdG = 0;
    double dIh2dG = 0;
    for (int i=0; i<I_list.size(); i++)
      {
      if (I_list[i].frame->ID == m)
        {
        d_numerator = I_list[i].value * I_list[i].sigma;
        d_denominator = 2 * I_list[i].sigma * I_list[i].frame->G;
        d_2_numerator = 0;
        d_2_denominator = 2 * I_list[i].sigma;
        }
      }
    dIhdG = (denominator * d_numerator - numerator * d_denominator) / (denominator * denominator);
    dIh2dG = ((denominator * denominator)*(d_numerator*d_denominator - (numerator * d_2_denominator+d_denominator*d_numerator))-(2 * denominator * d_denominator * (denominator * d_numerator - numerator * d_denominator))) / ((denominator * denominator)*(denominator * denominator));
    return derivatives(dIhdG, dIh2dG);

    }

};

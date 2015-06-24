#include <iostream>
#include <vector>
#include "scaling_data_structure.h"
using namespace std;
typedef scitbx::af::shared<int> shared_int;
typedef scitbx::af::shared<size_t> shared_size_t;
typedef scitbx::af::shared<double> shared_double;
typedef scitbx::af::shared<Miller> shared_miller;
typedef scitbx::af::shared<Frame> shared_frame;
typedef scitbx::af::versa<double, scitbx::af::flex_grid<> >  flex_grid_double;
int
kdelta(int m, int n) {
  if(m==n) return 1;
  return 0;
}

//////////////get the first and second drivatives///////////
static boost::python::tuple
get_f_g(shared_miller& MI, const shared_double& x, shared_frame& Fr) {
/////////////////////////////////////////////////
  cout<< "Assigning the refined paramters to Eqn. 1 "<<'\n';
  for (int i=0; i<Fr.size(); i++) {
       Fr[i].G = x[i];
    }

  cout<<"Done"<<'\n';
////////////////////////////////////////////////
  double f = 0;
  double t = 0;
  double Ih=0;
  cout<< "Computing the Target Function....."<<'\n';
  for (int i=0; i<MI.size(); i++)  {
    Ih=MI[i].Ih();
    for (int j=0;j<MI[i].I_list.size();j++)    {
      t = MI[i].I_list[j].value - (MI[i].I_list[j].frame->G * Ih);
      f = f+ MI[i].I_list[j].sigma *t * t;
    }
  }
  cout<< "Done.....The new value for F is: "<<f<<'\n';
////////////////////////////////////////////////
  cout<< "Calculating the first and second derivatives with respect to all scaling factors (Gm's)....."<<'\n';
  shared_double g(x.size(),0);
  shared_double c(x.size(),0);
  double dtdG = 0;
  double dt2dG2 = 0;
  double dIhdG=0;
  double dh2dG2=0;
  int delta=0;
  derivatives derv(0,0);
  for(int m=0; m<Fr.size();m++)  {   //Loop over all frames
    for (int h=0; h<Fr[m].miller_list.size(); h++) { //Loop over the miller indices that have dependencies on the m frame (i.e. have observations in that frame)
      Ih=Fr[m].miller_list[h]->Im;
      derv=Fr[m].miller_list[h]->dIhdG(m);
      dIhdG=derv.first_dirv;
      dh2dG2=derv.second_dirv;
  //    cout<<"First Derv "<<dIhdG<<" "<<"Second Derv "<<dh2dG2<<endl; // test by Nick to make sure no zeros
      for (int j=0;j<Fr[m].miller_list[h]->I_list.size();j++) { //Loop over the observations of the h miller index
        t = Fr[m].miller_list[h]->I_list[j].value - (Fr[m].miller_list[h]->I_list[j].frame->G * Ih);   //calculating the derivatives
        delta=kdelta(m, Fr[m].miller_list[h]->I_list[j].frame->ID);                                    // please see the latex document
        dtdG = -1 * ((Fr[m].miller_list[h]->I_list[j].frame->G * dIhdG) + (Ih*delta));
        g[m] += (2 * t * dtdG * Fr[m].miller_list[h]->I_list[j].sigma);
        dt2dG2 = -1 * ((Fr[m].miller_list[h]->I_list[j].frame->G * dh2dG2) + (2 * dIhdG*delta));
        c[m] +=   2 * (t * dt2dG2 + dtdG * dtdG) * Fr[m].miller_list[h]->I_list[j].sigma;
      }
    }
  }
  cout<<"Done"<<'\n';
  return (boost::python::make_tuple(f, g, c));
}

//////////////////read the intensties in miller object//////////

static boost::python::tuple
read_objects(shared_miller& MI) {
  shared_double Ih(MI.size(),0);
  shared_int M(MI.size(),0);
  for (int i=0; i<MI.size(); i++) {
    Miller m=MI[i];
//    cout<< "Miller ID "<< m.index<<" "<<m.Im<<'\n';
    Ih[i]=m.Im;
    M[i]=m.index;
  }

  return   (boost::python::make_tuple(Ih, M));
}



static boost::python::tuple 
calculate_residual_jacobian_wI_G(const shared_double& I,
                const shared_double& w,
                const shared_size_t& hkl,
                const shared_size_t& frames,
                const shared_double& S,
                const shared_double& G,
                const shared_double& Ih,
                const int n_frames,
                flex_grid_double& jacobian) {
 double Gm=0;
 double f=0;
 int m=0;
 int h=0;
 double dtdIh=0;
 double t = 0;
 shared_double residuals(I.size(),0);
    for (int i=0;i<I.size();i++)
    {
     m=frames[i];
     h = hkl[i];
     residuals[i] = w[i]*(I[i] - (G[m] * Ih[h] * S[h]));
     jacobian(i, m) = -w[i] * Ih[h] * S[h];
     jacobian(i, h+n_frames) = -w[i] * G[m] * S[h];
    }
  return (boost::python::make_tuple(residuals, jacobian));
}


///////////////////////////////This function Loads the miller and frame objects//////////////////////
static boost::python::object
Load_object(shared_double I,
            shared_double sigma,
            shared_int hkl,
            shared_int allMiller,
            shared_int frames, int nframes) {

  shared_frame Fr(nframes); //array that carries frame objects of all frames
  for(int i=0; i<nframes; i++) { //make a frame object to each frame in your data set
    Frame f;
    f.G=1;   //initialize the scale factor G to 1
    f.ID=i;
    Fr[i]=f;// add the new frame to the list of frames
  }
  int index=0;
  std::vector< Miller > MI;
  for (int i=0; i<allMiller.size(); i++) {
     Miller m;
     m.index=allMiller[i];
     for (int j=index; j<hkl.size(); j++) {
        if(hkl[j]==allMiller[i]) {
          Intensity In;
          In.value=I[j];
          In.sigma=sigma[j];
          In.frame= & Fr[frames[j]]; //assign the address of the frame object to the frame pointer based on the frame ID
          m.add_I(In);
          index++;
        }
        else
          break;
     }
       if (m.I_list.size()>0)
       MI.push_back(m);
  }

  shared_miller M(MI.size());
  for (int i=0; i<MI.size(); i++) {
  M[i]=MI[i];
  for (int j=0;j<M[i].I_list.size(); j++)
      if(!Fr[M[i].I_list[j].frame->ID].exist(& M[i]))
        Fr[M[i].I_list[j].frame->ID].miller_list.push_back(& M[i]);
  }
  return  (boost::python::make_tuple(M, Fr));
}

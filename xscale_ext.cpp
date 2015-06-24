#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <vector>
#include <map>
#include "minimizers/scaling_functions.cpp"
using namespace boost::python;
namespace xfel{
namespace boost_python { namespace {

  void
  xscale_module() {
    using namespace boost::python;
    def("Load_object", &Load_object);
    def("read_objects", &read_objects);
    def("get_f_g", &get_f_g);
    def("calculate_residual_jacobian_wI_G", &calculate_residual_jacobian_wI_G);
    typedef return_internal_reference<> rir;
    scitbx::af::boost_python::shared_wrapper<Miller, rir>::wrap("Miller");
    scitbx::af::boost_python::shared_wrapper<Frame, rir>::wrap("Frame");

  }

}
}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xscale_ext)
{
  xfel::boost_python::xscale_module();

}

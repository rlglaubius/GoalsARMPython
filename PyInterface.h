#ifndef PY_INTERFACE_H
#define PY_INTERFACE_H

#include <GoalsARM_Core.H>

#define BOOST_PYTHON_STATIC_LIB
#define BOOST_NUMPY_STATIC_LIB

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace py = boost::python;
namespace np = boost::python::numpy;

class PyInterface {
public:
	PyInterface(const int year_start, const int year_final);
	~PyInterface();

	inline void initialize(const std::string& upd_filename);
	void init_pasfrs_from_5yr(np::ndarray& pasfrs5y);
private:
	DP::Projection* proj;
};

BOOST_PYTHON_MODULE(GoalsARM) {
	np::initialize();

	py::class_<PyInterface>("Projection", py::init<size_t, size_t>())
		.def("initialize", &PyInterface::initialize)
		.def("init_pasfrs_from_5yr", &PyInterface::init_pasfrs_from_5yr)
	;
}

#endif // PY_INTERFACE_H

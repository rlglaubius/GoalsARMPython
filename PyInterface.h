#ifndef PY_INTERFACE_H
#define PY_INTERFACE_H

#include <GoalsARM_Core.H>

#define BOOST_PYTHON_STATIC_LIB
#define BOOST_NUMPY_STATIC_LIB

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace py = boost::python;
namespace np = boost::python::numpy;

// The PyInterface may take numpy ndarrays (np::ndarray) as input.
// However, PyInterface should never store pointers to any part of
// those arrays beyond the lifetime of the method that takes them
// as arguments. Python owns those arrays. If PyInterface does not
// handle reference counting for the arrays correctly, Python
// may garbage collect them, leaving PyInterface with dangling pointers
// that will eventually cause an error.

// TODO: We generally use boost::multi_array_ref to index into the data
// underlying np::ndarray. Since np::ndarray could be C-ordered or Fortran-ordered,
// these multi_array_ref instances are dynamically allocated. We should
// add exception handling to these methods that triggers when the
// ndarray size or arrangement is not compatible with the array consumer.
// Currently, these methods may cause an ugly crash instead.
// 
// See: https://docs.python.org/3/extending/extending.html

class PyInterface {
public:
	PyInterface(const int year_start, const int year_final);
	~PyInterface();

	inline void initialize(const std::string& upd_filename);

	void init_pasfrs_from_5yr(np::ndarray& pasfrs5y);
	void init_migr_from_5yr(np::ndarray& netmigr, np::ndarray& pattern_female, np::ndarray& pattern_male);

	// Toggle use of direct incidence. If flag=TRUE, use init_direct_incidence
	// to initialize incidence
	void use_direct_incidence(const bool flag);

	// initialize direct incidence inputs
	void init_direct_incidence(
		np::ndarray& inci,           // 1970-2050
		np::ndarray& sex_irr,        // 1970-2050
		np::ndarray& age_irr_female, // 1970-2050 by 5y age group (17 groups: 0-4, 5-9, ..., 75-79, 80+)
		np::ndarray& age_irr_male,   // 1970-2050 by 5y age group
		np::ndarray& pop_irr_female, // 1970-2050 by DP::N_POP
		np::ndarray& pop_irr_male);  // 1970-2050 by DP::N_POP

	void init_median_age_debut(const double age_female, const double age_male);
	void init_median_age_union(const double age_female, const double age_male);
	void init_mean_duration_union(const double years);

	// Turnover inputs
	void init_mean_duration_pwid(const double years_female, const double years_male);
	void init_mean_duration_fsw(const double years);
	void init_mean_duration_msm(const double years);

	// inputs are the proportion of people who are members of each key population
	// (e.g., the proportion of 15-49 females who are PWID, the proportion of 15+
	// males [assigned at birth] who are transgender women)
	void init_size_pwid(const double prop_female, const double prop_male);
	void init_size_fsw(const double prop);
	void init_size_msm(const double prop);
	void init_size_trans(const double prop_female, const double prop_male);

	// Initialize the age distribution of key populations with turnover
	void init_age_pwid(const double loc_female, const double shp_female, const double loc_male, const double shp_male);
	void init_age_fsw(const double loc, const double shp);
	void init_age_msm(const double loc, const double shp);

	// Initialize adult HIV progression and mortality rates off ART
	void init_adult_prog_from_10yr(np::ndarray& dist, np::ndarray& prog, np::ndarray& mort);

	// Initialize adult HIV-related mortality rates on ART
	void init_adult_art_mort_from_10yr(np::ndarray& art1, np::ndarray& art2, np::ndarray& art3, np::ndarray& art_mrr);


private:
	DP::Projection* proj;
};

BOOST_PYTHON_MODULE(GoalsARM) {
	np::initialize();

	py::class_<PyInterface>("Projection", py::init<size_t, size_t>())
		.def("initialize",               &PyInterface::initialize)
		.def("init_pasfrs_from_5yr",     &PyInterface::init_pasfrs_from_5yr)
		.def("init_migr_from_5yr",       &PyInterface::init_migr_from_5yr)
		.def("init_direct_incidence",    &PyInterface::init_direct_incidence)
		.def("init_median_age_debut",    &PyInterface::init_median_age_debut)
		.def("init_median_age_union",    &PyInterface::init_median_age_union)
		.def("init_mean_duration_union", &PyInterface::init_mean_duration_union)
		.def("init_mean_duration_pwid",  &PyInterface::init_mean_duration_pwid)
		.def("init_mean_duration_fsw",   &PyInterface::init_mean_duration_fsw)
		.def("init_mean_duration_msm",   &PyInterface::init_mean_duration_msm)
		.def("init_size_pwid",           &PyInterface::init_size_pwid)
		.def("init_size_fsw",            &PyInterface::init_size_fsw)
		.def("init_size_msm",            &PyInterface::init_size_msm)
		.def("init_size_trans",          &PyInterface::init_size_trans)
		.def("init_age_pwid",            &PyInterface::init_age_pwid)
		.def("init_age_fsw",             &PyInterface::init_age_fsw)
		.def("init_age_msm",             &PyInterface::init_age_msm)
		.def("init_adult_prog_from_10yr",     &PyInterface::init_adult_prog_from_10yr)
		.def("init_adult_art_mort_from_10yr", &PyInterface::init_adult_art_mort_from_10yr)

		.def("use_direct_incidence",     &PyInterface::use_direct_incidence)
	;
}

#endif // PY_INTERFACE_H

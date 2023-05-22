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

	// Initialize the first year of epidemic simulation, and HIV prevalence in that year
	void init_epidemic_seed(const int seed_year, const double seed_prev);

	// Initialize transmission probabilities per sex act
	void init_transmission(const double pct_f2m,
						   const double or_m2f,
						   const double or_m2m,
		                   const double primary,
						   const double chronic,
						   const double symptom,
						   const double art_supp,
						   const double art_fail);

	void init_hiv_fertility(np::ndarray& frr_age_off_art, np::ndarray& frr_cd4_off_art, np::ndarray& frr_age_on_art);

	// Initialize adult HIV progression and mortality rates off ART
	void init_adult_prog_from_10yr(np::ndarray& dist, np::ndarray& prog, np::ndarray& mort);

	// Initialize adult HIV-related mortality rates on ART
	void init_adult_art_mort_from_10yr(np::ndarray& art1, np::ndarray& art2, np::ndarray& art3, np::ndarray& art_mrr);

	void init_adult_art_eligibility(np::ndarray& cd4);
	void init_adult_art_curr(np::ndarray& art_num, np::ndarray& art_pct);

	// Initialize the ART initiation weight. ART uptake in a CD4
	// category is a weighted average of (1) the number of PLHIV off ART
	// in that category and (2) the expected number of deaths in that
	// population. "weight" is the weight assigned to expected deaths,
	// expressed as a percentage between 0 and 100.
	void init_adult_art_allocation(const double weight);
	
	// Initialize annual adult ART dropout rates from annual percentages
	// If 5% drop out each year, then the dropout rate is -ln(1-0.05)
	void init_adult_art_dropout(np::ndarray& art_drop_pct);

	// Initialize trends in adult viral suppression on ART
	// art_supp_pct is expected to have one row per year, and 8 columns. The
	// first four colums are for males ages 15-24, 25-34, 35-44, then 45+. The
	// next four columns are for females in those same age groups
	void init_adult_art_suppressed(np::ndarray& art_supp_pct);

	// Initialize male circumcision uptake
	// uptake is expected to have one row per year and one column per 5-year age group
	// 0-4, 5-9, ..., 75-79, 80+
	void init_male_circumcision_uptake(np::ndarray& uptake);

	void init_effect_vmmc(const double effect);
	void init_effect_condom(const double effect);

	// Initialize 14-year-old CLHIV from direct inputs
	void init_clhiv_agein(np::ndarray& clhiv);

	// Calculate projection from year_start to year_final
	void project(const int year_final);

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
		.def("init_epidemic_seed",       &PyInterface::init_epidemic_seed)
		.def("init_hiv_fertility",       &PyInterface::init_hiv_fertility)
		.def("init_transmission",        &PyInterface::init_transmission)
		.def("init_adult_prog_from_10yr",     &PyInterface::init_adult_prog_from_10yr)
		.def("init_adult_art_mort_from_10yr", &PyInterface::init_adult_art_mort_from_10yr)
		.def("init_adult_art_eligibility",    &PyInterface::init_adult_art_eligibility)
		.def("init_adult_art_curr",           &PyInterface::init_adult_art_curr)
		.def("init_adult_art_allocation",     &PyInterface::init_adult_art_allocation)
		.def("init_adult_art_dropout",        &PyInterface::init_adult_art_dropout)
		.def("init_adult_art_suppressed",     &PyInterface::init_adult_art_suppressed)
		.def("init_male_circumcision_uptake", &PyInterface::init_male_circumcision_uptake)
		.def("init_clhiv_agein",   &PyInterface::init_clhiv_agein)
		.def("init_effect_vmmc",   &PyInterface::init_effect_vmmc)
		.def("init_effect_condom", &PyInterface::init_effect_condom)

		.def("project", &PyInterface::project)

		.def("use_direct_incidence",     &PyInterface::use_direct_incidence)
	;
}

#endif // PY_INTERFACE_H

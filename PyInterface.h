#ifndef PY_INTERFACE_H
#define PY_INTERFACE_H

#include <GoalsARM.H>

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

	/// Pass memory for storing output population sizes.
	/// @param adult_neg HIV-negative adults, by year, sex, age, risk
	/// @param adult_hiv HIV-positive adults, by year, sex, age, risk, CD4, and care status
	/// @param child_neg HIV-negative children, by year, sex, age
	/// @param child_hiv HIV-positive children, by year, sex, age, CD4, and care status
	void share_output_population(
			np::ndarray& adult_neg,
			np::ndarray& adult_hiv,
			np::ndarray& child_neg,
			np::ndarray& child_hiv);

	/// Pass memory for storing output births counts.
	/// @param births Births by year and sex
	void share_output_births(np::ndarray& births);

	/// Pass memory for storing output births to mothers living with HIV
	/// @param births Births by year
	void share_output_births_exposed(np::ndarray& births);

	/// Pass memory for storing output all-cause deaths counts.
	/// @param adult_neg HIV-negative adults, by year, sex, age (15:80), risk
	/// @param adult_hiv HIV-positive adults, by year, sex, age (15:80), risk, CD4, and care status
	/// @param child_neg HIV-negative children, by year, sex, age (0:14)
	/// @param child_hiv HIV-positive children, by year, sex, age (0:14), CD4, and care status
	/// @details sex should have three levels: females, uncircumcised males, circumcised males
	void share_output_deaths(
			np::ndarray& adult_neg,
			np::ndarray& adult_hiv,
			np::ndarray& child_neg,
			np::ndarray& child_hiv);

	/// Pass memory for storing output new HIV infections
	/// @param newhiv New HIV infections by year, sex, age (0:80), risk
	/// @details sex should have three levels: females, uncircumcised males, circumcised males
	void share_output_new_infections(np::ndarray& newhiv);

	/// Pass partner rate inputs
	/// @param partner_rate matrix by year (year_start:year_final), sex (male,female), age (15:80), and behavioral risk group
	void share_input_partner_rate(np::ndarray& partner_rate);

	/// Pass mixing preferences by age
	/// @param mix Mixing matrix. mix[si,ai,sj,aj] is the preference of sex si, age ai for partners of sex sj, age aj
	/// Generally, mix[si,ai,sj,:] should sum to 1. Age indices 0..66 correspond to ages 15..80. Mixing coefficients
	/// for female-female partnerships or age 80 are not currently used.
	void share_input_age_mixing(np::ndarray& mix);

	/// Pass assortativity parameters for behavioral risk groups
	/// @param assort assort[s,r] is the extent that people of sex s and behavioral risk r mix preferentially
	void share_input_pop_assort(np::ndarray& assort);

	/// Pass parameters that specify HIV acquisition risk in people who inject drugs
	/// @param force Force of infection acting on PWID who share needles, by year and sex
	/// @param needle_sharing proportion of PWID who share needles by year
	void share_input_pwid_risk(np::ndarray& force, np::ndarray& needle_sharing);

	/// Use a UPD file to initialize demographic inputs
	/// @param upd_filename UPD file name
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

	/// Initialize key population size parameter values
	/// @param kp_size kp_size[i] is the proportion of the overall the 15-49 population in key population i.
	/// @param kp_stay kp_stay[i] is 1 if people stay in the key population after entry, 0 if they may eventually leave.
	/// @param kp_turnover kp_turnover[k,i] indices the average number of years spent in the population (k=0),
	/// the median age of population members (k=1) and the age distribution shape parameter (k=2). These may be left
	/// uninitialized if people remain in the population after entry.
	/// @details
	/// Populations:
	/// - i=0 people who inject drugs, female
	/// - i=1 people who inject drugs, male
	/// - i=2 female sex workers
	/// - i=3 male clients of female sex workers
	/// - i=4 men who have sex with men
	/// - i=5 transgender women
	void init_keypop_size_params(
		np::ndarray& kp_size,
		np::ndarray& kp_stay,
		np::ndarray& kp_turnover);

	/// Initialize the proportion of key population members who have a main opposite-sex partner
	/// @param prop_married prop_married[i] is the proportion of key population i who are married to or
	/// have a cohabiting opposite-sex partner.
	/// @details
	/// Populations:
	/// - i=0 people who inject drugs, female
	/// - i=1 people who inject drugs, male
	/// - i=2 female sex workers
	/// - i=3 male clients of female sex workers
	/// - i=4 men who have sex with men
	/// - i=5 transgender women
	void init_keypop_married(np::ndarray& prop_married);

	/// Initialize the structure of the mixing matrix by behavioral risk group
	/// @param mix_levels matrix of mixing levels. mix_levels[si,ri,sj,rj] for
	/// (sex, risk group) pairs (si,ri) and (sj,rj) takes values 0, 1, or 2. These
	/// indicate if the groups do not mix (value=0), can mix (value=1) or prefer
	/// to mix (value=2). We do not require that this matrix is symmetric. Sexes
	/// si and sj refer to assigned sex at birth, not to gender identity.
	void init_mixing_matrix(np::ndarray& mix_levels);

	/// Initialize numbers of sex acts per year by partnership type
	/// @param acts a vector storing the number of sex acts per year by
	/// partnership type (main=0, casual=1, commercial=2, msm=3)
	void init_sex_acts(np::ndarray& acts);
	
	/// Initialize condom use inputs by year and partnership type
	/// @param freq freq[t][i] is the probability in [0,1] of condom use at last sex
	/// by partnership type i (main=0, casual=1, commercial=2, msm=3)
	void init_condom_freq(np::ndarray& freq);

	/// Initialize the first year of epidemic simulation, and HIV prevalence in that year
	/// @param seed_year First year of the HIV epidemic. This should be specified as the number of years since the projection began
	/// @param seed_prev HIV prevalence in the first year of the HIV epidemic.
	void init_epidemic_seed(const int seed_year, const double seed_prev);

	/// Initialize transmission probabilities per sex act
	/// @param transmit_f2m female-to-male transmission probability (as proportion) per sex act
	/// @param or_m2f odds ratio for male-to-female transmission, relative to female-to-male
	/// @param or_m2m odds ratio for male-to-male transmission, relative to female-to-male
	/// @param primary odds ratio for transmission during primary infection
	/// @param chronic odds ratio for transmission during chronic (asymtomatic) infection
	/// @param symptom odds ratio for transmission during symptomatic infection
	/// @param or_art_supp odds ratio for transmission on ART when virally suppressed, relative to off ART
	/// @param or_art_fail odds ratio for transmission on ART when virally unsuppressed, relative to off ART
	void init_transmission(
			const double transmit_f2m,
			const double or_m2f,
			const double or_m2m,
			const double primary,
			const double chronic,
			const double symptom,
			const double or_art_supp,
			const double or_art_fail);

	void init_hiv_fertility(np::ndarray& frr_age_off_art, np::ndarray& frr_cd4_off_art, np::ndarray& frr_age_on_art);

	// Initialize adult HIV progression and mortality rates off ART
	void init_adult_prog_from_10yr(np::ndarray& dist, np::ndarray& prog, np::ndarray& mort);

	// Initialize adult HIV-related mortality rates on ART
	void init_adult_art_mort_from_10yr(np::ndarray& art1, np::ndarray& art2, np::ndarray& art3, np::ndarray& art_mrr);

	/// Initialize CD4 thresholds for adult ART eligibility
	/// @param cd4 CD4 thresholds over time. This must be of integer type.
	void init_adult_art_eligibility(np::ndarray& cd4);

	/// Initialize adult ART program size
	/// @param n_art number on ART by year and sex (female, male)
	/// @param p_art proportion in [0,1] of ART need met by year and sex (female, male)
	/// @details Adult ART coverage can be specified in absolute numbers or as a proportion
	/// of need met. For a given year t and sex s, if the proportion p_art[t][s] > 0 then it
	/// is used to drive calculations, otherwise n_art[t][s] is used.
	void init_adult_art_curr(np::ndarray& n_art, np::ndarray& p_art);

	// Initialize the ART initiation weight. ART uptake in a CD4
	// category is a weighted average of (1) the number of PLHIV off ART
	// in that category and (2) the expected number of deaths in that
	// population. "weight" is the weight assigned to expected deaths,
	// expressed as a percentage between 0 and 100.
	void init_adult_art_allocation(const double weight);
	
	/// Initialize annual adult ART interruption rates
	/// @param art_exit_rate ART interruption rates by year
	/// @details This should be an event rate (interruptions per person-year),
	/// not a proportion or percentage
	void init_adult_art_interruption(np::ndarray& art_exit_rate);

	// Initialize trends in adult viral suppression on ART
	// art_supp_pct is expected to have one row per year, and 8 columns. The
	// first four colums are for males ages 15-24, 25-34, 35-44, then 45+. The
	// next four columns are for females in those same age groups
	void init_adult_art_suppressed(np::ndarray& art_supp_pct);

	// Initialize male circumcision uptake
	// uptake is expected to have one row per year and one column per 5-year age group
	// 0-4, 5-9, ..., 75-79, 80+
	void init_male_circumcision_uptake(np::ndarray& uptake);

	/// Initialize the effect of male circumcision on HIV acquisition
	/// @param effect proportionate reduction (in [0,1]) in HIV acquisition when circumcised (vs. not)
	/// @details This is used for direct and for mechanistic incidence calculations
	void init_effect_vmmc(const double effect);

	/// Initialize the effect of condom use on HIV transmission
	/// @param effect proportionate reduction (in [0,1]) in HIV transmission per act when a condom is used (vs. not)
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
		.def("share_output_population",     &PyInterface::share_output_population)
		.def("share_output_births",         &PyInterface::share_output_births)
		.def("share_output_deaths",         &PyInterface::share_output_deaths)
		.def("share_output_births_exposed", &PyInterface::share_output_births_exposed)
		.def("share_output_new_infections", &PyInterface::share_output_new_infections)
		.def("share_input_partner_rate",    &PyInterface::share_input_partner_rate)
		.def("share_input_age_mixing",      &PyInterface::share_input_age_mixing)
		.def("share_input_pop_assort",      &PyInterface::share_input_pop_assort)
		.def("share_input_pwid_risk",       &PyInterface::share_input_pwid_risk)

		.def("initialize",               &PyInterface::initialize)
		.def("init_pasfrs_from_5yr",     &PyInterface::init_pasfrs_from_5yr)
		.def("init_migr_from_5yr",       &PyInterface::init_migr_from_5yr)
		.def("init_direct_incidence",    &PyInterface::init_direct_incidence)
		.def("init_median_age_debut",    &PyInterface::init_median_age_debut)
		.def("init_median_age_union",    &PyInterface::init_median_age_union)
		.def("init_mean_duration_union", &PyInterface::init_mean_duration_union)
		.def("init_keypop_size_params",  &PyInterface::init_keypop_size_params)
		.def("init_keypop_married",      &PyInterface::init_keypop_married)
		.def("init_mixing_matrix",       &PyInterface::init_mixing_matrix)
		.def("init_sex_acts",            &PyInterface::init_sex_acts)
		.def("init_condom_freq",         &PyInterface::init_condom_freq)
		.def("init_epidemic_seed",       &PyInterface::init_epidemic_seed)
		.def("init_hiv_fertility",       &PyInterface::init_hiv_fertility)
		.def("init_transmission",        &PyInterface::init_transmission)
		.def("init_adult_prog_from_10yr",     &PyInterface::init_adult_prog_from_10yr)
		.def("init_adult_art_mort_from_10yr", &PyInterface::init_adult_art_mort_from_10yr)
		.def("init_adult_art_eligibility",    &PyInterface::init_adult_art_eligibility)
		.def("init_adult_art_curr",           &PyInterface::init_adult_art_curr)
		.def("init_adult_art_allocation",     &PyInterface::init_adult_art_allocation)
		.def("init_adult_art_interruption",   &PyInterface::init_adult_art_interruption)
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

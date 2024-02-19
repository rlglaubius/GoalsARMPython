#ifndef GOALS_PROJ_H
#define GOALS_PROJ_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <GoalsARM.h>

namespace py = pybind11;

typedef py::array_t<int> array_int_t;
typedef py::array_t<double> array_double_t;

class GoalsProj {
public:
	GoalsProj(const int year_start, const int year_final);
	~GoalsProj();

	/// Pass memory for storing output population sizes.
	/// @param adult_neg HIV-negative adults, by year, sex, age, risk
	/// @param adult_hiv HIV-positive adults, by year, sex, age, risk, CD4, and care status
	/// @param child_neg HIV-negative children, by year, sex, age
	/// @param child_hiv HIV-positive children, by year, sex, age, CD4, and care status
	void share_output_population(
		array_double_t adult_neg,
		array_double_t adult_hiv,
		array_double_t child_neg,
		array_double_t child_hiv);

	/// Pass memory for storing output births counts.
	/// @param births Births by year and sex
	void share_output_births(array_double_t births);

	/// Pass memory for storing output all-cause deaths counts.
	/// @param adult_neg HIV-negative adults, by year, sex, age (15:80), risk
	/// @param adult_hiv HIV-positive adults, by year, sex, age (15:80), risk, CD4, and care status
	/// @param child_neg HIV-negative children, by year, sex, age (0:14)
	/// @param child_hiv HIV-positive children, by year, sex, age (0:14), CD4, and care status
	/// @details sex should have three levels: females, uncircumcised males, circumcised males
	void share_output_deaths(
		array_double_t adult_neg,
		array_double_t adult_hiv,
		array_double_t child_neg,
		array_double_t child_hiv);

	/// Pass memory for storing output new HIV infections
	/// @param newhiv New HIV infections by year, sex, age (0:80), risk
	/// @details sex should have three levels: females, uncircumcised males, circumcised males
	void share_output_new_infections(array_double_t newhiv);

	/// Pass memory for storing output births to mothers living with HIV
	/// @param births Births by year
	void share_output_births_exposed(array_double_t births);

	/// Pass partner rate inputs
	/// @param partner_rate matrix by year (year_start:year_final), sex (male,female), age (15:80), and behavioral risk group
	void share_input_partner_rate(array_double_t partner_rate);

	/// Pass mixing preferences by age
	/// @param mix Mixing matrix. mix[si,ai,sj,aj] is the preference of sex si, age ai for partners of sex sj, age aj
	/// Generally, mix[si,ai,sj,:] should sum to 1. Age indices 0..66 correspond to ages 15..80. Mixing coefficients
	/// for female-female partnerships or age 80 are not currently used.
	void share_input_age_mixing(array_double_t mix);

	/// Pass assortativity parameters for behavioral risk groups
	/// @param assort Array by sex and behavioral risk group. assort[s,r] is
	/// the extent that people of sex s and behavioral risk r mix preferentially
	void share_input_pop_assort(array_double_t assort);

	/// Pass parameters that specify HIV acquisition risk in people who inject drugs
	/// @param force Force of infection acting on PWID who share needles, by year and sex
	/// @param needle_sharing proportion of PWID who share needles by year
	void share_input_pwid_risk(array_double_t force, array_double_t needle_sharing);

	/// Use a UPD file to initialize demographic inputs
	/// @param upd_filename UPD file name
	void initialize(const std::string& upd_filename);

	/// Initialize proportionate age-specific fertility (PASFR) from inputs by five-year age group
	/// @param pasfrs5y an array by year and age group (15-19, 20-24, ..., 45-49)
	/// This initialization method is provided for compatibility with Spectrum.
	void init_pasfrs_from_5yr(array_double_t pasfrs5y);

	/// Initialize net migration from inputs by five-year age group
	/// @param netmigr total  net migrants by year and sex
	/// @param pattern_female proportionate migration of females by year and age group (0-4, 5-9, ..., 75-79, 80+)
	/// @param pattern_male   proportionate migration of males by year and age group (0-4, 5-9, ..., 75-79, 80+)
	/// Proportionate migration is calculated as the absolute number of net migrants in an age group, divided by the
	/// overall number of net migrants. This initialization method is provided for compatibility with Spectrum.
	void init_migr_from_5yr(array_double_t netmigr, array_double_t pattern_female, array_double_t pattern_male);

	/// Initialize direct incidence inputs
	/// @param inci      array of incidence rates (infections per person-year) by year
	/// @param sex_irr   array of incidence rate ratios for females relative to males by year
	/// @param age_irr_f array of incidence rate ratios for females by year and five-year age group (0-4, 5-9, ..., 75-79, 80+)
	/// @param age_irr_m array of incidence rate ratios for males by year and five-year age group
	/// @param pop_irr_f array of incidence rate ratios for females by behavioral risk group
	/// @param pop_irr_m array of incidence rate ratios for males by behavioral risk group
	void init_direct_incidence(
		array_double_t inci,
		array_double_t sex_irr,
		array_double_t age_irr_f,
		array_double_t age_irr_m,
		array_double_t pop_irr_f,
		array_double_t pop_irr_m);

	/// Initialize the median age at sexual debut
	/// @param age_female median age at sexual debut among females
	/// @param age_male   median age at sexual debut among males
	void init_median_age_debut(const double age_female, const double age_male);

	/// Initialize the median age at first union (marriage or cohabitation)
	/// @param age_female median age at first union among females
	/// @param age_male   median age at first union among males
	void init_median_age_union(const double age_female, const double age_male);

	/// Initialize the average duration of marriage or cohabitation in years
	/// @param years average duration
	void init_mean_duration_union(const double years);

	/// Initialize key population size parameter values
	/// @param kp_size kp_size[i] is the proportion of the overall the 15-49 population in key population i.
	/// @param kp_stay kp_stay[i] is 1 if people stay in the key population after entry, 0 if they may eventually leave.
	/// @param kp_turnover kp_turnover[k,i] stores the average number of years spent in the population (k=0),
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
	void init_keypop_size_params(array_double_t kp_size, array_int_t kp_stay, array_double_t kp_turnover);

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
	void init_keypop_married(array_double_t prop_married);

	/// Initialize the structure of the mixing matrix by behavioral risk group
	/// @param mix_levels matrix of mixing levels. mix_levels[si,ri,sj,rj] for
	/// (sex, risk group) pairs (si,ri) and (sj,rj) takes values 0, 1, or 2. These
	/// indicate if the groups do not mix (value=0), can mix (value=1) or prefer
	/// to mix (value=2). We do not require that this matrix is symmetric. Sexes
	/// si and sj refer to assigned sex at birth, not to gender identity.
	void init_mixing_matrix(array_double_t mix_levels);

	/// Initialize numbers of sex acts per year by partnership type
	/// @param acts a vector storing the number of sex acts per year by
	/// partnership type (main=0, casual=1, commercial=2, msm=3)
	void init_sex_acts(array_double_t acts);

	/// Initialize condom use inputs by year and partnership type
	/// @param freq freq[t][i] is the probability in [0,1] of condom use at last sex
	/// by partnership type i (main=0, casual=1, commercial=2, msm=3)
	void init_condom_freq(array_double_t freq);

	/// Initialize input STI symptom prevalence trends
	/// @param sti_prev Array by year, sex, age, and behavioral risk group
	void init_sti_prev(array_double_t sti_prev);

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
	/// @param or_sti_hiv_pos odds ratio for HIV transmission for STI symptoms in HIV-positive partner
	/// @param or_sti_hiv_neg odds ratio for HIV transmission for STI symptoms in HIV-negative partner
	void init_transmission(
		const double transmit_f2m,
		const double or_m2f,
		const double or_m2m,
		const double primary,
		const double chronic,
		const double symptom,
		const double or_art_supp,
		const double or_art_fail,
		const double or_sti_hiv_pos,
		const double or_sti_hiv_neg);

	/// Initialize HIV-related fertility rate ratios (FRRs)
	/// @param frr_age_off_art FRRs off ART by year and five-year age group (15-19, 20-24, ..., 45-49)
	/// @param frr_cd4_off_art FRRs off ART by HIV infection stage (primary, CD4>500, CD4 350-500, ..., CD4<50)
	/// @param frr_age_on_art  FRRs on ART by five-year age group
	void init_hiv_fertility(array_double_t frr_age_off_art, array_double_t frr_cd4_off_art, array_double_t frr_age_on_art);

	/// Initialize adult HIV progression and mortality rates off ART
	/// @param dist HIV stage at infection by CD4 category, excluding primary infection
	/// @param prog HIV disease progression rates by HIV stage, excluding the last, since prog=0 is implied
	/// @param mort HIV mortality rates by HIV stage
	/// @details dist, prog and mort must be 2-d arrays with 8 columns corresponding to sex, age combinations:
	/// male 15-24, male 25-34, male 35-44, male 45+, female 15-24, female 25-34, female 35-44, female 45+
	void init_adult_prog_from_10yr(array_double_t dist, array_double_t prog, array_double_t mort);

	/// Initialize adult HIV-related mortality rates on ART
	/// @param art1 HIV-related mortality rates by CD4 category when on ART for [0,6) months
	/// @param art2 HIV-related mortality rates by CD4 category when on ART for [6,12) months
	/// @param art3 HIV-related mortality rates by CD4 category when on ART for 12+ months
	/// @param art_mrr HIV-related mortality rate ratios on ART by year and time on ART ([0,12), 12+ months)
	/// @details art1, art2 and art3 must be 2-d arrays with 8 columns corresponding to sex, age combinations:
	/// male 15-24, male 25-34, male 35-44, male 45+, female 15-24, female 25-34, female 35-44, female 45+
	void init_adult_art_mort_from_10yr(array_double_t art1, array_double_t art2, array_double_t art3, array_double_t art_mrr);

	/// Initialize CD4 thresholds for adult ART eligibility
	/// @param cd4 CD4 thresholds by year
	void init_adult_art_eligibility(array_int_t cd4);

	/// Initialize adult ART program size
	/// @param n_art number on ART by year and sex (female, male)
	/// @param p_art proportion in [0,1] of ART need met by year and sex (female, male)
	/// @details Adult ART coverage can be specified in absolute numbers or as a proportion
	/// of need met. For a given year t and sex s, if the proportion p_art[t][s] > 0 then it
	/// is used to drive calculations, otherwise n_art[t][s] is used.
	void init_adult_art_curr(array_double_t n_art, array_double_t p_art);

	/// Initialize the ART initiation weight.
	/// @param weight ART initiation weight.
	/// @details ART uptake in a CD4 category is a weighted average of (1) the number of
	/// PLHIV off ART in that category and (2) the expected number of deaths in that
	/// population. "weight" is the weight assigned to expected deaths, expressed as a
	/// proportion between 0 and 1.
	void init_adult_art_allocation(const double weight);

	/// Initialize annual adult ART interruption rates
	/// @param art_exit_rate ART interruption rates by year and sex
	/// @details This should be an event rate (interruptions per person-year),
	/// not a proportion or percentage
	void init_adult_art_interruption(array_double_t art_exit_rate);

	/// Initialize trends in adult viral suppression on ART
	/// @param art_supp_prop proportion of adults on ART who are virally suppressed by year, age, and sex
	/// @details art_supp_prop must have 8 columns corresponding to sex, age combinations:
	/// male 15-24, male 25-34, male 35-44, male 45+, female 15-24, female 25-34, female 35-44, female 45+
	void init_adult_art_suppressed(array_double_t art_supp_pct);

	/// Initialize male circumcision uptake
	/// @param uptake uptake rates by year and five-year age group (0-4, 5-9, ..., 75-79, 80+)
	void init_male_circumcision_uptake(array_double_t uptake);

	/// Initialize the effect of male circumcision on HIV acquisition
	/// @param effect proportionate reduction (in [0,1]) in HIV acquisition when circumcised (vs. not)
	/// @details This is used for direct and for mechanistic incidence calculations
	void init_effect_vmmc(const double effect);

	/// Initialize the effect of condom use on HIV transmission
	/// @param effect proportionate reduction (in [0,1]) in HIV transmission per act when a condom is used (vs. not)
	void init_effect_condom(const double effect);

	/// Initialize 14-year-old CLHIV from direct inputs
	/// @param clhiv a 2-d array of children living with HIV by year
	/// @details columns of CLHIV correspond to Spectrum strata: sex, pediatric CD4 category,
	/// ART duration ([0,6), [6,12), 12+ months) if on ART or HIV acquisition timing
	/// (perinatal, breastfeeding within [0,6), [6,12), 12+ months of birth). Since Spectrum
	/// forgets transmission timing once children start ART, there should be 84 rows (2 sexes,
	/// 6 CD4, 7 ART durations or acquisition timings).
	void init_clhiv_agein(array_double_t clhiv);

	/// Calculate the projection
	/// @param year_final the last year to project
	/// @details If project(...) is called repeatedly, each calculation will
	/// resume from the latest year calculated in previous calls. Use invalidate(...)
	/// to resume calculations from an earlier year.
	void project(const int year_final);


	/// Invalidate projected calculations from year onward.
	/// @param year Invalidate calculated outcomes from this year onward.
	/// @details After project(t) is called, subsequent calls to project(...)
	/// will not recalculate years <= t. Use invalidate(...) to reset 
	/// this to a selected year. Setting year < 0 will cause the next
	/// project(...) call to start from the first year of projection.
	inline void invalidate(const int year);

	/// Toggle use of direct incidence. If flag=TRUE, 
	/// @param flag If TRUE if direct incidence inputs should be used,
	/// FALSE if mechanistic incidence calculations should be done
	void use_direct_incidence(const bool flag);

private:
	DP::Projection* proj;
	size_t num_years;
};

// GoalsProj is an interface to the calculation engine
// 
// Memory management: The calculation engine uses some workspaces allocated by
// Python client applications. The client retains ownership of this memory.
// Methods for passing variables to the calculation engine fall into three
// categories: share_input, share_output, and init.
// 
// share_input: The calculation engine retains references to arguments for the
// lifetime of the GoalsProj instance. These arguments are not modified by the
// calculation engine.
// 
// share_output: The calculation engine retains references to arguments for the
// lifetime of the GoalsProj instance. These arguments may be modified by the
// calculation engine.
// 
// init: The calculation engine does not modify or retain references to arguments.
// Arguments are safe to garbage collect after the method returns control to the
// client.
// 
// py::keep_alive<1,n>() keeps the garbage collector from deallocating argument 
// n so long as the GoalsProj instance (argument 1) is still alive
PYBIND11_MODULE(goals_proj, m) {
	py::class_<GoalsProj>(m, "Projection")
		.def(py::init<const int, const int>())

		.def("share_output_population",     &GoalsProj::share_output_population,     py::keep_alive<1,2>(), py::keep_alive<1,3>(), py::keep_alive<1,4>(), py::keep_alive<1,5>())
		.def("share_output_births",         &GoalsProj::share_output_births,         py::keep_alive<1,2>())
		.def("share_output_deaths",         &GoalsProj::share_output_deaths,         py::keep_alive<1,2>(), py::keep_alive<1,3>(), py::keep_alive<1,4>(), py::keep_alive<1,5>())
		.def("share_output_births_exposed", &GoalsProj::share_output_births_exposed, py::keep_alive<1,2>())
		.def("share_output_new_infections", &GoalsProj::share_output_new_infections, py::keep_alive<1,2>())
		.def("share_input_partner_rate",    &GoalsProj::share_input_partner_rate,    py::keep_alive<1,2>())
		.def("share_input_age_mixing",      &GoalsProj::share_input_age_mixing,      py::keep_alive<1,2>())
		.def("share_input_pop_assort",	    &GoalsProj::share_input_pop_assort,      py::keep_alive<1,2>())
		.def("share_input_pwid_risk",       &GoalsProj::share_input_pwid_risk,       py::keep_alive<1,2>(), py::keep_alive<1,3>())

		.def("initialize",                    &GoalsProj::initialize)
		.def("init_pasfrs_from_5yr",          &GoalsProj::init_pasfrs_from_5yr)
		.def("init_migr_from_5yr",            &GoalsProj::init_migr_from_5yr)
		.def("init_direct_incidence",         &GoalsProj::init_direct_incidence)
		.def("init_median_age_debut",         &GoalsProj::init_median_age_debut)
		.def("init_median_age_union",         &GoalsProj::init_median_age_union)
		.def("init_mean_duration_union",      &GoalsProj::init_mean_duration_union)
		.def("init_keypop_size_params",       &GoalsProj::init_keypop_size_params)
		.def("init_keypop_married",           &GoalsProj::init_keypop_married)
		.def("init_mixing_matrix",            &GoalsProj::init_mixing_matrix)
		.def("init_sex_acts",                 &GoalsProj::init_sex_acts)
		.def("init_condom_freq",              &GoalsProj::init_condom_freq)
		.def("init_sti_prev",                 &GoalsProj::init_sti_prev)
		.def("init_epidemic_seed",            &GoalsProj::init_epidemic_seed)
		.def("init_hiv_fertility",            &GoalsProj::init_hiv_fertility)
		.def("init_transmission",             &GoalsProj::init_transmission)
		.def("init_adult_prog_from_10yr",     &GoalsProj::init_adult_prog_from_10yr)
		.def("init_adult_art_mort_from_10yr", &GoalsProj::init_adult_art_mort_from_10yr)
		.def("init_adult_art_eligibility",    &GoalsProj::init_adult_art_eligibility)
		.def("init_adult_art_curr",           &GoalsProj::init_adult_art_curr)
		.def("init_adult_art_allocation",     &GoalsProj::init_adult_art_allocation)
		.def("init_adult_art_interruption",   &GoalsProj::init_adult_art_interruption)
		.def("init_adult_art_suppressed",     &GoalsProj::init_adult_art_suppressed)
		.def("init_male_circumcision_uptake", &GoalsProj::init_male_circumcision_uptake)
		.def("init_clhiv_agein",              &GoalsProj::init_clhiv_agein)
		.def("init_effect_vmmc",              &GoalsProj::init_effect_vmmc)
		.def("init_effect_condom",            &GoalsProj::init_effect_condom)

		.def("project",    &GoalsProj::project)
		.def("invalidate", &GoalsProj::invalidate)

		.def("use_direct_incidence", &GoalsProj::use_direct_incidence)

		;

#ifdef VERSION_INFO
#else
	m.attr("__version__") = "dev";
#endif
}

#endif // GOALS_PROJ_H

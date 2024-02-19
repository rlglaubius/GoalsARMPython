#include <format>
#include <boost/math/interpolators/pchip.hpp>
#include "goals_proj.h"

template<typename ValueType>
ValueType* prepare_array(py::array_t<ValueType> arr, const size_t ndim, size_t* shape) {
	auto buff = arr.request();

	if (buff.ndim != ndim) {
		throw std::runtime_error(std::format("array dimension={}, expected {}", buff.ndim, ndim));
	}

	for (int dim(0); dim < ndim; ++dim) {
		if (buff.shape[dim] != shape[dim]) {
			throw std::runtime_error(std::format("array dimension {} has {} elements, expected {}", dim, buff.shape[dim], shape[dim]));
		}
	}

	if (!(arr.flags() & py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_)) {
		throw std::runtime_error(std::format("array must be a contiguous C-style array"));
	}

	if (!(arr.flags() & py::detail::npy_api::constants::NPY_ARRAY_ALIGNED_)) {
		throw std::runtime_error(std::format("array is not properly aligned"));
	}

	return reinterpret_cast<ValueType*>(buff.ptr);
}

GoalsProj::GoalsProj(const int year_start, const int year_final)
	: num_years(year_final - year_start + 1) {
	proj = new DP::Projection(year_start, year_final);
}

GoalsProj::~GoalsProj() {
	if (proj != NULL) { delete proj; }
}

void GoalsProj::share_output_population(
	array_double_t adult_neg,
	array_double_t adult_hiv,
	array_double_t child_neg,
	array_double_t child_hiv) {
	size_t shape_adult_neg[] = {num_years, DP::N_SEX_MC, DP::N_AGE_ADULT, DP::N_POP};
	size_t shape_adult_hiv[] = {num_years, DP::N_SEX_MC, DP::N_AGE_ADULT, DP::N_POP, DP::N_HIV_ADULT, DP::N_DTX};
	size_t shape_child_neg[] = {num_years, DP::N_SEX_MC, DP::N_AGE_CHILD};
	size_t shape_child_hiv[] = {num_years, DP::N_SEX_MC, DP::N_AGE_CHILD, DP::N_HIV_CHILD, DP::N_DTX};
	proj->pop.share_storage(
		prepare_array(adult_neg, 4, shape_adult_neg),
		prepare_array(adult_hiv, 6, shape_adult_hiv),
		prepare_array(child_neg, 3, shape_child_neg),
		prepare_array(child_hiv, 5, shape_child_hiv));
}

void GoalsProj::share_output_births(array_double_t births) {
	size_t shape[] = {num_years, DP::N_SEX};
	double* ptr_births(prepare_array(births, 2, shape));
	proj->dat.share_births(ptr_births);
}

void GoalsProj::share_output_deaths(
	array_double_t adult_neg,
	array_double_t adult_hiv,
	array_double_t child_neg,
	array_double_t child_hiv) {
	size_t shape_adult_neg[] = {num_years, DP::N_SEX_MC, DP::N_AGE_ADULT, DP::N_POP};
	size_t shape_adult_hiv[] = {num_years, DP::N_SEX_MC, DP::N_AGE_ADULT, DP::N_POP, DP::N_HIV_ADULT, DP::N_DTX};
	size_t shape_child_neg[] = {num_years, DP::N_SEX_MC, DP::N_AGE_CHILD};
	size_t shape_child_hiv[] = {num_years, DP::N_SEX_MC, DP::N_AGE_CHILD, DP::N_HIV_CHILD, DP::N_DTX};
	proj->dth.share_storage(
		prepare_array(adult_neg, 4, shape_adult_neg),
		prepare_array(adult_hiv, 6, shape_adult_hiv),
		prepare_array(child_neg, 3, shape_child_neg),
		prepare_array(child_hiv, 5, shape_child_hiv));
}

void GoalsProj::share_output_new_infections(array_double_t newhiv) {
	size_t shape[] = {num_years, DP::N_SEX_MC, DP::N_AGE, DP::N_POP};
	double* ptr_newhiv(prepare_array(newhiv, 4, shape));
	proj->dat.share_new_infections(ptr_newhiv);
}

void GoalsProj::share_output_births_exposed(array_double_t births) {
	size_t shape[] = {num_years};
	double* ptr_births(prepare_array(births, 1, shape));
	proj->dat.share_births_exposed(ptr_births);
}

void GoalsProj::share_input_partner_rate(array_double_t partner_rate) {
	size_t shape[] = {num_years, DP::N_SEX, DP::N_AGE_ADULT, DP::N_POP};
	double* ptr_partner_rate(prepare_array(partner_rate, 4, shape));
	proj->dat.share_partner_rate(ptr_partner_rate);
}

void GoalsProj::share_input_age_mixing(array_double_t mix) {
	size_t shape[] = {DP::N_SEX, DP::N_AGE_ADULT, DP::N_SEX, DP::N_AGE_ADULT};
	double* ptr_mix(prepare_array(mix, 4, shape));
	proj->dat.share_age_mixing(ptr_mix);
}

void GoalsProj::share_input_pop_assort(array_double_t assort) {
	size_t shape[] = {DP::N_SEX, DP::N_POP};
	double* ptr_assort(prepare_array(assort, 2, shape));
	proj->dat.share_pop_assortativity(ptr_assort);
}

void GoalsProj::share_input_pwid_risk(array_double_t force, array_double_t needle_sharing) {
	size_t shape_force[] = {num_years, DP::N_SEX};
	size_t shape_share[] = {num_years};
	double* ptr_force(prepare_array(force, 2, shape_force));
	double* ptr_share(prepare_array(needle_sharing, 1, shape_share));
	proj->dat.share_pwid_risk(ptr_force, ptr_share);
}

void GoalsProj::initialize(const std::string& upd_filename) {
	proj->initialize(upd_filename);
}

void GoalsProj::init_pasfrs_from_5yr(array_double_t pasfrs5y) {
	size_t shape[] = {num_years, 7};
	double* ptr_fert(prepare_array(pasfrs5y, 2, shape));
	DP::year_age_ref_t fert(ptr_fert, boost::extents[shape[0]][shape[1]]);
	proj->dat.init_pasfrs_from_5yr(fert);
}

void GoalsProj::init_migr_from_5yr(array_double_t netmigr, array_double_t pattern_female, array_double_t pattern_male) {
	const int n_age5y(17);
	size_t shape_netmigr[] = {num_years, DP::N_SEX};
	size_t shape_pattern[] = {num_years, n_age5y};

	double* ptr_netmigr(prepare_array(netmigr, 2, shape_netmigr));
	double* ptr_pattern_f(prepare_array(pattern_female, 2, shape_pattern));
	double* ptr_pattern_m(prepare_array(pattern_male,   2, shape_pattern));

	// We use boost containers here just to simplify indexing into the flat arrays
	// underlying the input arrays.
	DP::year_sex_ref_t migr(ptr_netmigr, boost::extents[shape_netmigr[0]][shape_netmigr[1]]);
	DP::year_age_ref_t migr_f(ptr_pattern_f, boost::extents[shape_pattern[0]][shape_pattern[1]]);
	DP::year_age_ref_t migr_m(ptr_pattern_m, boost::extents[shape_pattern[0]][shape_pattern[1]]);

	// Convert patterns from multipliers to absolute net migrant numbers
	for (int t(0); t < num_years; ++t) {
		for (int a(0); a < n_age5y; ++a) {
			migr_m[t][a] = migr[t][0] * migr_m[t][a]; // Sex order is swapped in Excel (Male, Female) compared to Goals ARM (Female, Male)
			migr_f[t][a] = migr[t][1] * migr_f[t][a];
		}
	}

	proj->dat.init_migr_from_5yr(DP::FEMALE, migr_f);
	proj->dat.init_migr_from_5yr(DP::MALE,   migr_m);
}

void GoalsProj::init_direct_incidence(
	array_double_t inci,
	array_double_t sex_irr,
	array_double_t age_irr_f,
	array_double_t age_irr_m,
	array_double_t pop_irr_f,
	array_double_t pop_irr_m) {

	const int n_age5y(17);
	size_t shape_age[] = {num_years, n_age5y};
	size_t shape_pop[] = {num_years, DP::N_POP};
	size_t shape_series[] = {num_years};

	double* ptr_inci(prepare_array(inci, 1, shape_series));
	double* ptr_sex_irr(prepare_array(sex_irr, 1, shape_series));
	double* ptr_age_irr_f(prepare_array(age_irr_f, 2, shape_age));
	double* ptr_age_irr_m(prepare_array(age_irr_m, 2, shape_age));
	double* ptr_pop_irr_f(prepare_array(pop_irr_f, 2, shape_pop));
	double* ptr_pop_irr_m(prepare_array(pop_irr_m, 2, shape_pop));

	DP::year_age_ref_t airr_f(ptr_age_irr_f, boost::extents[shape_age[0]][shape_age[1]]);
	DP::year_age_ref_t airr_m(ptr_age_irr_m, boost::extents[shape_age[0]][shape_age[1]]);
	boost::multi_array_ref<double, 2> pirr_f(ptr_pop_irr_f, boost::extents[shape_pop[0]][shape_pop[1]]);
	boost::multi_array_ref<double, 2> pirr_m(ptr_pop_irr_m, boost::extents[shape_pop[0]][shape_pop[1]]);

	for (int t(0); t < num_years; ++t) proj->dat.incidence(t, ptr_inci[t]);
	for (int t(0); t < num_years; ++t) proj->dat.irr_sex(t, ptr_sex_irr[t]);
	for (int t(0); t < num_years; ++t) {
		for (int r(DP::POP_MIN); r <= DP::POP_MAX; ++r) {
			proj->dat.irr_pop(t, DP::FEMALE, r, pirr_f[t][r]);
			proj->dat.irr_pop(t, DP::MALE,   r, pirr_m[t][r]);
		}
	}

	proj->dat.init_age_irr_from_5yr(DP::FEMALE, airr_f);
	proj->dat.init_age_irr_from_5yr(DP::MALE,   airr_m);
}

void GoalsProj::init_median_age_debut(const double age_female, const double age_male) {
	DP::set_median_age_debut(proj->dat, DP::FEMALE, age_female);
	DP::set_median_age_debut(proj->dat, DP::MALE,   age_male);
}

void GoalsProj::init_median_age_union(const double age_female, const double age_male) {
	DP::set_median_age_union(proj->dat, DP::FEMALE, age_female);
	DP::set_median_age_union(proj->dat, DP::MALE,   age_male);
}

void GoalsProj::init_mean_duration_union(const double years) {
	DP::set_mean_union_duration(proj->dat, years);
}

void GoalsProj::init_keypop_size_params(array_double_t kp_size, array_int_t kp_stay, array_double_t kp_turnover) {
	const int n_pop(6), n_elt(3);
	const DP::sex_t sex[] = {DP::FEMALE,   DP::MALE,     DP::FEMALE,  DP::MALE,    DP::MALE,    DP::MALE};
	const DP::pop_t pop[] = {DP::POP_PWID, DP::POP_PWID, DP::POP_FSW, DP::POP_CSW, DP::POP_MSM, DP::POP_TGW};
	bool stay;
	double med, loc, shp;

	size_t shape[] = {n_pop};
	size_t shape_turn[] = {n_elt, n_pop};
	int* ptr_stay(prepare_array(kp_stay, 1, shape));
	double* ptr_size(prepare_array(kp_size, 1, shape));
	double* ptr_turn(prepare_array(kp_turnover, 2, shape_turn));
	boost::multi_array_ref<double,2> turnover(ptr_turn, boost::extents[shape_turn[0]][shape_turn[1]]);

	for (int r(0); r < n_pop; ++r) {
		stay = ptr_stay[r];
		proj->dat.keypop_size(sex[r], pop[r], ptr_size[r]);
		proj->dat.keypop_stay(sex[r], pop[r], stay);
		if (!stay) {
			med = turnover[1][r];
			shp = turnover[2][r];
			loc = log(med - 15.0);
			DP::set_mean_keypop_duration(proj->dat, sex[r], pop[r], turnover[0][r]);
			DP::set_keypop_age(proj->dat, sex[r], pop[r], loc, shp);
		}
	}
}

void GoalsProj::init_keypop_married(array_double_t prop_married) {
	const int n_pop(6);
	const DP::sex_t sex[] = {DP::FEMALE,   DP::MALE,     DP::FEMALE,  DP::MALE,    DP::MALE,    DP::MALE};
	const DP::pop_t pop[] = {DP::POP_PWID, DP::POP_PWID, DP::POP_FSW, DP::POP_CSW, DP::POP_MSM, DP::POP_TGW};
	size_t shape[] = {n_pop};
	double* ptr_married(prepare_array(prop_married, 1, shape));
	for (int r(0); r < n_pop; ++r) {
		proj->dat.keypop_married(sex[r], pop[r], ptr_married[r]);
	}
}

void GoalsProj::init_mixing_matrix(array_double_t mix_levels) {
	size_t shape[] = {DP::N_SEX, DP::N_POP, DP::N_SEX, DP::N_POP};
	double* ptr_mix_levels(prepare_array(mix_levels, 4, shape));
	boost::multi_array_ref<double, 4> arr_mix_levels(ptr_mix_levels, boost::extents[shape[0]][shape[1]][shape[2]][shape[3]]);
	for (int si(DP::SEX_MIN); si <= DP::SEX_MAX; ++si)
		for (int ri(DP::POP_MIN); ri <= DP::POP_MAX; ++ri)
			for (int sj(DP::SEX_MIN); sj <= DP::SEX_MAX; ++sj)
				for (int rj(DP::POP_MIN); rj <= DP::POP_MAX; ++rj)
					proj->dat.mix_structure(si, ri, sj, rj, arr_mix_levels[si][ri][sj][rj]);
}

void GoalsProj::init_sex_acts(array_double_t acts) {
	size_t shape[] = {DP::N_BOND};
	double* ptr_acts(prepare_array(acts, 1, shape));
	for (int q(DP::BOND_MIN); q <= DP::BOND_MAX; ++q)
		proj->dat.sex_acts(q, ptr_acts[q]);
}

void GoalsProj::init_condom_freq(array_double_t freq) {
	size_t shape[] = {num_years, DP::N_BOND};
	double* ptr_freq(prepare_array(freq, 2, shape));
	boost::multi_array_ref<double, 2> arr_freq(ptr_freq, boost::extents[shape[0]][shape[1]]);
	for (int t(0); t < proj->num_years(); ++t)
		for (int q(DP::BOND_MIN); q <= DP::BOND_MAX; ++q)
			proj->dat.condom_freq(t, q, arr_freq[t][q]);
}

void GoalsProj::init_sti_prev(array_double_t sti_prev) {
	const int ndim(4);
	size_t shape[] = {num_years, DP::N_SEX, DP::N_AGE_ADULT, DP::N_POP};
	double* ptr_sti_prev(prepare_array(sti_prev, ndim, shape));
	boost::multi_array_ref<double, ndim> arr_sti_prev(ptr_sti_prev, boost::extents[shape[0]][shape[1]][shape[2]][shape[3]]);
	for (int t(0); t < shape[0]; ++t)
		for (int s(0); s < shape[1]; ++s)
			for (int a(0); a < shape[2]; ++a)
				for (int r(0); r < shape[3]; ++r)
					proj->dat.sti_prev(t, s, a, r, arr_sti_prev[t][s][a][r]);
}

void GoalsProj::init_epidemic_seed(const int seed_year, const double seed_prev) {
	proj->dat.seed_time(seed_year);
	proj->dat.seed_prevalence(seed_prev);
}

void GoalsProj::init_hiv_fertility(array_double_t frr_age_off_art, array_double_t frr_cd4_off_art, array_double_t frr_age_on_art) {
	const int n_age5y(7); // 15-19, 20-24, ..., 45-49
	size_t shape_age[] = {num_years, n_age5y};
	size_t shape_cd4[] = {DP::N_HIV};
	size_t shape_art[] = {n_age5y};

	double* ptr_age(prepare_array(frr_age_off_art, 2, shape_age));
	double* frr_cd4(prepare_array(frr_cd4_off_art, 1, shape_cd4));
	double* frr_art(prepare_array(frr_age_on_art,  1, shape_art));
	boost::multi_array_ref frr_age(ptr_age, boost::extents[shape_age[0]][shape_age[1]]);

	int b;
	for (int a(0); a < DP::N_AGE_BIRTH; ++a) {
		b = a / 5;
		for (int t(0); t < proj->num_years(); ++t)
			proj->dat.frr_age_no_art(t, a, frr_age[t][b]);
		proj->dat.frr_age_on_art(a, frr_art[b]);
	}
	for (int h(DP::HIV_ADULT_MIN); h <= DP::HIV_ADULT_MAX; ++h)
		proj->dat.frr_cd4_no_art(h, frr_cd4[h]);
}

void GoalsProj::init_transmission(
	const double transmit_f2m,
	const double or_m2f,
	const double or_m2m,
	const double primary,
	const double chronic,
	const double symptom,
	const double or_art_supp,
	const double or_art_fail,
	const double or_sti_hiv_pos,
	const double or_sti_hiv_neg) {
	DP::set_transmission(proj->dat,
		transmit_f2m,
		or_m2f,
		or_m2m,
		primary,
		chronic,
		symptom,
		or_art_supp,
		or_art_fail,
		or_sti_hiv_pos,
		or_sti_hiv_neg);
}

void GoalsProj::init_adult_prog_from_10yr(array_double_t dist, array_double_t prog, array_double_t mort) {
	size_t shape_prog[] = {DP::N_HIV - 1, DP::N_SEX * 4};
	size_t shape_mort[] = {DP::N_HIV,     DP::N_SEX * 4};
	DP::cd4_sex_age_ref_t arr_dist(prepare_array(dist, 2, shape_prog), boost::extents[shape_prog[0]][shape_prog[1]]);
	DP::cd4_sex_age_ref_t arr_prog(prepare_array(prog, 2, shape_prog), boost::extents[shape_prog[0]][shape_prog[1]]);
	DP::cd4_sex_age_ref_t arr_mort(prepare_array(mort, 2, shape_mort), boost::extents[shape_mort[0]][shape_mort[1]]);

	DP::set_adult_prog_from_10yr(proj->dat, arr_dist, arr_prog, arr_mort);
}

void GoalsProj::init_adult_art_mort_from_10yr(array_double_t art1, array_double_t art2, array_double_t art3, array_double_t art_mrr) {
	size_t shape_art[] = {DP::N_HIV, DP::N_SEX * 4};
	size_t shape_mrr[] = {num_years, 2}; // This uses 2 ART durations ([0,12), [12,\infty) months on ART) instead of Goals's 3 ([0,6), [6,12), [12,\infty)).

	DP::cd4_sex_age_ref_t arr_art1(prepare_array(art1, 2, shape_art), boost::extents[shape_art[0]][shape_art[1]]);
	DP::cd4_sex_age_ref_t arr_art2(prepare_array(art2, 2, shape_art), boost::extents[shape_art[0]][shape_art[1]]);
	DP::cd4_sex_age_ref_t arr_art3(prepare_array(art3, 2, shape_art), boost::extents[shape_art[0]][shape_art[1]]);
	DP::year_dtx_ref_t arr_mrr(prepare_array(art_mrr, 2, shape_mrr), boost::extents[shape_mrr[0]][shape_mrr[1]]);

	DP::set_adult_art_mort_from_10yr(proj->dat, arr_art1, arr_art2, arr_art3, arr_mrr);
}

void GoalsProj::init_adult_art_eligibility(array_int_t cd4) {
	size_t shape[] = {num_years};
	DP::time_series_int_ref_t arr_cd4(prepare_array(cd4, 1, shape), boost::extents[shape[0]]);
	DP::set_adult_art_eligibility_from_cd4(proj->dat, arr_cd4);
}

void GoalsProj::init_adult_art_curr(array_double_t n_art, array_double_t p_art) {
	size_t shape[] = {num_years, DP::N_SEX};
	DP::year_sex_ref_t arr_n_art(prepare_array(n_art, 2, shape), boost::extents[shape[0]][shape[1]]);
	DP::year_sex_ref_t arr_p_art(prepare_array(p_art, 2, shape), boost::extents[shape[0]][shape[1]]);
	for (int t(0); t < proj->dat.num_years(); ++t) {
		proj->dat.art_num_adult(t, DP::MALE,   arr_n_art[t][0]);
		proj->dat.art_num_adult(t, DP::FEMALE, arr_n_art[t][1]);
		proj->dat.art_prop_adult(t, DP::MALE,   arr_p_art[t][0]);
		proj->dat.art_prop_adult(t, DP::FEMALE, arr_p_art[t][1]);
	}
}

void GoalsProj::init_adult_art_allocation(const double weight) {
	proj->dat.art_mort_weight(weight);
}

void GoalsProj::init_adult_art_interruption(array_double_t art_exit_rate) {
	size_t shape[] = {num_years, DP::N_SEX};
	DP::year_sex_ref_t arr_exit(prepare_array(art_exit_rate, 2, shape), boost::extents[shape[0]][shape[1]]);
	for (int t(0); t < proj->dat.num_years(); ++t) {
		proj->dat.art_exit_adult(t, DP::MALE,   arr_exit[t][0]);
		proj->dat.art_exit_adult(t, DP::FEMALE, arr_exit[t][1]);
	}
}

void GoalsProj::init_adult_art_suppressed(array_double_t art_supp_pct) {
	const int n_age(4);
	size_t shape[] = {num_years, DP::N_SEX * n_age};
	boost::multi_array_ref<double, 2> arr_supp(prepare_array(art_supp_pct, 2, shape), boost::extents[shape[0]][shape[1]]);
	int col_m, col_f;
	for (int t(0); t < proj->dat.num_years(); ++t) {
		for (int a(0); a < DP::N_AGE_ADULT; ++a) {
			col_m = std::min(a / 10, n_age - 1);
			col_f = col_m + n_age;
			proj->dat.art_suppressed_adult(t, DP::MALE,   a, arr_supp[t][col_m]);
			proj->dat.art_suppressed_adult(t, DP::FEMALE, a, arr_supp[t][col_f]);
		}
	}
}

void GoalsProj::init_male_circumcision_uptake(array_double_t uptake) {
	using boost::math::interpolators::pchip;

	const size_t n(17); // number of 5-year age groups
	std::vector<double> x{0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85};
	std::vector<double> y(n + 1);
	double rate, prop;
	double dy_bgn, dy_end; // derivatives at left and right boundaries

	size_t shape[] = {num_years, n};
	DP::year_age_ref_t arr_uptake(prepare_array(uptake, 2, shape), boost::extents[shape[0]][shape[1]]);

	y[0] = 0.0;
	for (int t(0); t < proj->dat.num_years(); ++t) {
		// Calculate cumulative exposure to circumcision uptake at the
		// boundaries of five-year age groups
		for (int a(0); a < n; ++a) {
			prop = 0.01 * arr_uptake[t][a];
			rate = -5.0 * log(1.0 - prop);
			y[a+1] = y[a] + rate;
		}

		// Boundary calculations below are designed to match signal::pchip in
		// R, which uses public Fortran code by one of the authors of PCHIP
		// (doi:10.1137/0717021). R produces different results than boost::pchip
		// does by default because R approximates boundary derivatives using a
		// quadratic formula while boost::pchip uses a linear formula.
		// 
		// We have not implemented the general formula for calculating
		// derivatives, since we know ages in x are spaced five years apart (the
		// hard-coded constants 0.2, 1.5, and -0.5 come from that optimization).
		dy_bgn = 0.2 * (1.5 * (y[1] - y[ 0 ]) - 0.5 * (y[ 2 ] - y[ 1 ]));
		dy_end = 0.2 * (1.5 * (y[n] - y[n-1]) - 0.5 * (y[n-2] - y[n-1]));
		if (dy_bgn < 0.0) dy_bgn = 0.0;
		if (dy_end < 0.0) dy_end = 0.0;

		// Interpolate cumulative exposure at singles ages using PCHIP. We
		// pass copies of vectors to pchip because it is allowed to modify
		// its inputs, including resizing them
		auto spline = pchip(std::vector<double>(x), std::vector<double>(y), dy_bgn, dy_end);

		// Calculate incremental uptake between consecutive ages and convert
		// back from rates to proportions
		for (int a(0); a < DP::N_AGE; ++a) {
			rate = spline(a + 1) - spline(a);
			prop = 1.0 - exp(-rate);
			proj->dat.uptake_male_circumcision(t, a, prop);
		}
	}
}

void GoalsProj::init_effect_vmmc(const double effect) {
	proj->dat.effect_vmmc(effect);
}

void GoalsProj::init_effect_condom(const double effect) {
	proj->dat.effect_condom(effect);
}

void GoalsProj::init_clhiv_agein(array_double_t clhiv) {
	size_t shape[] = {num_years, DP::N_SEX * (DP::HIV_CHILD_PED_MAX - DP::HIV_CHILD_PED_MIN + 1) * 7};
	boost::multi_array_ref<double, 2> arr_clhiv(prepare_array(clhiv, 2, shape), boost::extents[shape[0]][shape[1]]);
	DP::set_clhiv_agein(proj->dat, arr_clhiv);
}

void GoalsProj::project(const int year_final) {
	proj->project(year_final);
}

void GoalsProj::invalidate(const int year) {
	proj->invalidate(year);
}

void GoalsProj::use_direct_incidence(const bool flag) {
	proj->dat.direct_incidence(flag);
}

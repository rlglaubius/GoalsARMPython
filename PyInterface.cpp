#include "PyInterface.h"
#include <fstream>
#include <iostream>
#include <boost/math/interpolators/pchip.hpp>

template<typename ValueType>
ValueType* prepare_ndarray(np::ndarray& arr) {
#ifdef _DEBUG // _DEBUG is set when compiling in Debug mode
		if (arr.get_nd() > 1 && !(arr.get_flags() & (np::ndarray::C_CONTIGUOUS | np::ndarray::ALIGNED))) {
				// numpy multidimensional arrays must be C contiguous arrays, otherwise
				// Goals calculation code will access elements out-of-order. This will generally
				// produce incorrect results, so we abort execution instead.
				//
				// How might this happen? Suppose A is a numpy array. It may have the wrong
				// memory ordering if:
				// - If A was created in Fortran order (e.g., numpy.array(..., order="F"))
				// - If A was created C order, then transposed via A.transpose().
				// - Other operations, like numpyp.reshape(), might also change flags instead of memory ordering
				throw std::runtime_error("Array misaligned");
		}
#endif // _DEBUG
		return reinterpret_cast<ValueType*>(arr.get_data());
}

// Factory functions. Each returns a pointer to an n-dimensional boost::multi_array_ref
// that is a wrapper around the data underlying an input numpy ndarray.
// 
// CAVEATS:
// The template parameter NumType must be implicitly castable from the arr dtype.
// The caller owns the pointer returned, and is responsible for deallocating it.
template<typename NumType> boost::multi_array_ref<NumType, 1>* generate1d(np::ndarray& arr);
template<typename NumType> boost::multi_array_ref<NumType, 2>* generate2d(np::ndarray& arr);

// TODO: Handle potential exception cases
// sizeof(arr.dtype()) != sizeof(NumType)
// ptr = nullptr
// arr.get_nd() != target dimension
// review other arr.get_flags() for potential issues

template<typename NumType>
boost::multi_array_ref<NumType, 1>* generate1d(np::ndarray& arr) {
	const int ext = arr.shape(0);
	NumType* data_ptr = reinterpret_cast<NumType*>(arr.get_data());
	boost::multi_array_ref<NumType, 1>* ptr;
	ptr = new boost::multi_array_ref<NumType, 1>(data_ptr, boost::extents[ext]);
	return ptr;
}

template<typename NumType>
boost::multi_array_ref<NumType, 2>* generate2d(np::ndarray& arr) {
	const int ext0 = arr.shape(0);
	const int ext1 = arr.shape(1);
	NumType* data_ptr = reinterpret_cast<NumType*>(arr.get_data());
	boost::multi_array_ref<NumType, 2>* ptr;

	if (arr.get_flags() & np::ndarray::C_CONTIGUOUS) {
		ptr = new boost::multi_array_ref<NumType,2>(data_ptr, boost::extents[ext0][ext1], boost::c_storage_order());
	} else if (arr.get_flags() & np::ndarray::F_CONTIGUOUS) {
		ptr = new boost::multi_array_ref<NumType, 2>(data_ptr, boost::extents[ext0][ext1], boost::fortran_storage_order());
	} else {
		ptr = nullptr;
	}

	return ptr;
}

PyInterface::PyInterface(const int year_start, const int year_final) {
	proj = new DP::Projection(year_start, year_final);
}

PyInterface::~PyInterface() {
	if (proj != NULL) { delete proj; }
}

void PyInterface::share_output_population(
		np::ndarray& adult_neg,
		np::ndarray& adult_hiv,
		np::ndarray& child_neg,
		np::ndarray& child_hiv) {
		double* ptr_adult_neg(prepare_ndarray<double>(adult_neg));
		double* ptr_adult_hiv(prepare_ndarray<double>(adult_hiv));
		double* ptr_child_neg(prepare_ndarray<double>(child_neg));
		double* ptr_child_hiv(prepare_ndarray<double>(child_hiv));
		proj->pop.share_storage(ptr_adult_neg, ptr_adult_hiv, ptr_child_neg, ptr_child_hiv);
}

void PyInterface::share_output_births(np::ndarray& births) {
	double* ptr_births(prepare_ndarray<double>(births));
	proj->dat.share_births(ptr_births);
}

void PyInterface::share_output_deaths(
		np::ndarray& adult_neg,
		np::ndarray& adult_hiv,
		np::ndarray& child_neg,
		np::ndarray& child_hiv) {
		double* ptr_adult_neg(prepare_ndarray<double>(adult_neg));
		double* ptr_adult_hiv(prepare_ndarray<double>(adult_hiv));
		double* ptr_child_neg(prepare_ndarray<double>(child_neg));
		double* ptr_child_hiv(prepare_ndarray<double>(child_hiv));
		proj->dth.share_storage(ptr_adult_neg, ptr_adult_hiv, ptr_child_neg, ptr_child_hiv);
}

void PyInterface::share_output_new_infections(np::ndarray& newhiv) {
	double* ptr_newhiv(prepare_ndarray<double>(newhiv));
	proj->dat.share_new_infections(ptr_newhiv);
}

void PyInterface::share_input_partner_rate(np::ndarray& partner_rate) {
		double* ptr_partner_rate(prepare_ndarray<double>(partner_rate));
		proj->dat.share_partner_rate(ptr_partner_rate);
}

void PyInterface::share_input_age_mixing(np::ndarray& mix) {
	double* ptr_mix(prepare_ndarray<double>(mix));
	proj->dat.share_age_mixing(ptr_mix);
}

void PyInterface::share_input_pop_assort(np::ndarray& assort) {
	double* ptr_assort(prepare_ndarray<double>(assort));
	proj->dat.share_pop_assortativity(ptr_assort);
}

void PyInterface::share_input_pwid_risk(np::ndarray& force, np::ndarray& needle_sharing) {
	double* ptr_force(prepare_ndarray<double>(force));
	double* ptr_share(prepare_ndarray<double>(needle_sharing));
	proj->dat.share_pwid_risk(ptr_force, ptr_share);
}

void PyInterface::initialize(const std::string& upd_filename) {
	proj->initialize(upd_filename);
}

void PyInterface::init_pasfrs_from_5yr(np::ndarray& pasfrs5y) {
	DP::year_age_ref_t* dat = generate2d<double>(pasfrs5y);
	if (dat != NULL) {
		proj->dat.init_pasfrs_from_5yr(*dat);
		delete dat;
	}
}

void PyInterface::init_migr_from_5yr(np::ndarray& netmigr, np::ndarray& pattern_female, np::ndarray& pattern_male) {
	DP::year_sex_ref_t* migr = generate2d<double>(netmigr);
	DP::year_age_ref_t* migr_f = generate2d<double>(pattern_female);
	DP::year_age_ref_t* migr_m = generate2d<double>(pattern_male);

	const size_t nt(migr_f->shape()[0]), na(migr_f->shape()[1]);

	// convert patterns from multipliers to absolute net migrant numbers
	// since migr, migr_m, and migr_f are pointers, we do have to do some
	// ugly dereferencing
	for (int t(0); t < nt; ++t)
		for (int a(0); a < na; ++a) {
			(*migr_m)[t][a] = (*migr)[t][0] * (*migr_m)[t][a];
			(*migr_f)[t][a] = (*migr)[t][1] * (*migr_f)[t][a];
		}

	proj->dat.init_migr_from_5yr(DP::FEMALE, *migr_f);
	proj->dat.init_migr_from_5yr(DP::MALE,   *migr_m);

	delete migr;
	delete migr_f;
	delete migr_m;
}

void PyInterface::use_direct_incidence(const bool flag) {
	proj->dat.direct_incidence(flag);
}

void PyInterface::init_direct_incidence(
	np::ndarray& inci,
	np::ndarray& sex_irr,
	np::ndarray& age_irr_female,
	np::ndarray& age_irr_male,
	np::ndarray& pop_irr_female,
	np::ndarray& pop_irr_male) {

	DP::year_age_ref_t* airr_f = generate2d<double>(age_irr_female);
	DP::year_age_ref_t* airr_m = generate2d<double>(age_irr_male);

	const int nt(proj->year_final() - proj->year_first() + 1);

	for (int t(0); t < nt; ++t) proj->dat.incidence(t, py::extract<double>(inci[t]));
	for (int t(0); t < nt; ++t) proj->dat.irr_sex(t, py::extract<double>(sex_irr[t]));
	for (int t(0); t < nt; ++t)
		for (int r(DP::POP_MIN); r <= DP::POP_MAX; ++r) {
			proj->dat.irr_pop(t, DP::FEMALE, r, py::extract<double>(pop_irr_female[t][r]));
			proj->dat.irr_pop(t, DP::MALE,   r, py::extract<double>(pop_irr_male[t][r]));
		}

	proj->dat.init_age_irr_from_5yr(DP::FEMALE, *airr_f);
	proj->dat.init_age_irr_from_5yr(DP::MALE,   *airr_m);
}

void PyInterface::init_median_age_debut(const double age_female, const double age_male) {
	DP::set_median_age_debut(proj->dat, DP::FEMALE, age_female);
	DP::set_median_age_debut(proj->dat, DP::MALE,   age_male);
}

void PyInterface::init_median_age_union(const double age_female, const double age_male) {
	DP::set_median_age_union(proj->dat, DP::FEMALE, age_female);
	DP::set_median_age_union(proj->dat, DP::MALE, age_male);
}

void PyInterface::init_mean_duration_union(const double years) {
	DP::set_mean_union_duration(proj->dat, years);
}

void PyInterface::init_keypop_size_params(
	np::ndarray& kp_size,
	np::ndarray& kp_stay,
	np::ndarray& kp_turnover) {

	const int n_pop(6);
	const DP::sex_t sex[] = {DP::FEMALE,   DP::MALE,     DP::FEMALE,  DP::MALE,    DP::MALE,    DP::MALE};
	const DP::pop_t pop[] = {DP::POP_PWID, DP::POP_PWID, DP::POP_FSW, DP::POP_CSW, DP::POP_MSM, DP::POP_TGW};
	bool stay;
	double med, loc, shp;

	for (int r(0); r < n_pop; ++r) {
		stay = py::extract<int>(kp_stay[r]);
		proj->dat.keypop_size(sex[r], pop[r], py::extract<double>(kp_size[r]));
		proj->dat.keypop_stay(sex[r], pop[r], stay);
		if (!stay) {
			med = py::extract<double>(kp_turnover[1][r]);
			shp = py::extract<double>(kp_turnover[2][r]);
			loc = log(med - 15.0);
			DP::set_mean_keypop_duration(proj->dat, sex[r], pop[r], py::extract<double>(kp_turnover[0][r]));
			DP::set_keypop_age(proj->dat, sex[r], pop[r], loc, shp);
		}
	}
}

void PyInterface::init_keypop_married(np::ndarray& prop_married) {
	const int n_pop(6);
	const DP::sex_t sex[] = {DP::FEMALE,   DP::MALE,     DP::FEMALE,  DP::MALE,    DP::MALE,    DP::MALE};
	const DP::pop_t pop[] = {DP::POP_PWID, DP::POP_PWID, DP::POP_FSW, DP::POP_CSW, DP::POP_MSM, DP::POP_TGW};
	for (int r(0); r < n_pop; ++r) {
		proj->dat.keypop_married(sex[r], pop[r], py::extract<double>(prop_married[r]));
	}
}

void PyInterface::init_mixing_matrix(np::ndarray& mix_levels) {
	for (int si(DP::SEX_MIN); si <= DP::SEX_MAX; ++si)
		for (int ri(DP::POP_MIN); ri <= DP::POP_MAX; ++ri)
			for (int sj(DP::SEX_MIN); sj <= DP::SEX_MAX; ++sj)
				for (int rj(DP::POP_MIN); rj <= DP::POP_MAX; ++rj)
					proj->dat.mix_structure(si, ri, sj, rj, py::extract<int>(mix_levels[si][ri][sj][rj]));
}

void PyInterface::init_sex_acts(np::ndarray& acts) {
	for (int q(DP::BOND_MIN); q <= DP::BOND_MAX; ++q)
		proj->dat.sex_acts(q, py::extract<double>(acts[q]));
}

void PyInterface::init_condom_freq(np::ndarray& freq) {
	for (int t(0); t < proj->num_years(); ++t)
		for (int q(DP::BOND_MIN); q <= DP::BOND_MAX; ++q)
			proj->dat.condom_freq(t, q, py::extract<double>(freq[t][q]));
}

void PyInterface::init_epidemic_seed(const int seed_year, const double seed_prev) {
	proj->dat.seed_time(seed_year);
	proj->dat.seed_prevalence(seed_prev);
}

void PyInterface::init_hiv_fertility(np::ndarray& frr_age_off_art, np::ndarray& frr_cd4_off_art, np::ndarray& frr_age_on_art) {
	int b;
	for (int a(0); a < DP::N_AGE_BIRTH; ++a) {
		b = a / 5;
		for (int t(0); t < proj->num_years(); ++t)
			proj->dat.frr_age_no_art(t, a, py::extract<double>(frr_age_off_art[t][b]));
		proj->dat.frr_age_on_art(a, py::extract<double>(frr_age_on_art[b]));
	}
	for (int h(DP::HIV_ADULT_MIN); h <= DP::HIV_ADULT_MAX; ++h)
		proj->dat.frr_cd4_no_art(h, py::extract<double>(frr_cd4_off_art[h]));
}

void PyInterface::init_transmission(
	const double transmit_f2m,
	const double or_m2f,
	const double or_m2m,
	const double primary,
	const double chronic,
	const double symptom,
	const double or_art_supp,
	const double or_art_fail) {
	DP::set_transmission(proj->dat, transmit_f2m, or_m2f, or_m2m, primary, chronic, symptom, or_art_supp, or_art_fail);
}

void PyInterface::init_adult_prog_from_10yr(np::ndarray& dist, np::ndarray& prog, np::ndarray& mort) {
	DP::cd4_sex_age_ref_t* ptr_dist = generate2d<double>(dist);
	DP::cd4_sex_age_ref_t* ptr_prog = generate2d<double>(prog);
	DP::cd4_sex_age_ref_t* ptr_mort = generate2d<double>(mort);

	DP::set_adult_prog_from_10yr(proj->dat, *ptr_dist, *ptr_prog, *ptr_mort);

	delete ptr_dist;
	delete ptr_prog;
	delete ptr_mort;
}

void PyInterface::init_adult_art_mort_from_10yr(np::ndarray& art1, np::ndarray& art2, np::ndarray& art3, np::ndarray& art_mrr) {
	DP::cd4_sex_age_ref_t* ptr_art1 = generate2d<double>(art1);
	DP::cd4_sex_age_ref_t* ptr_art2 = generate2d<double>(art2);
	DP::cd4_sex_age_ref_t* ptr_art3 = generate2d<double>(art3);
	DP::year_dtx_ref_t* ptr_mrr = generate2d<double>(art_mrr);

	DP::set_adult_art_mort_from_10yr(proj->dat, *ptr_art1, *ptr_art2, *ptr_art3, *ptr_mrr);

	delete ptr_art1;
	delete ptr_art2;
	delete ptr_art3;
	delete ptr_mrr;
}

void PyInterface::init_adult_art_eligibility(np::ndarray& cd4) {
	DP::time_series_int_ref_t* ptr_cd4 = generate1d<int>(cd4);
	DP::set_adult_art_eligibility_from_cd4(proj->dat, *ptr_cd4);
	delete ptr_cd4;
}

void PyInterface::init_adult_art_curr(np::ndarray& n_art, np::ndarray& p_art) {
	for (int t(0); t < proj->dat.num_years(); ++t) {
		proj->dat.art_num_adult(t, DP::MALE,   py::extract<double>(n_art[t][0]));
		proj->dat.art_num_adult(t, DP::FEMALE, py::extract<double>(n_art[t][1]));
		proj->dat.art_prop_adult(t, DP::MALE,   py::extract<double>(p_art[t][0]));
		proj->dat.art_prop_adult(t, DP::FEMALE, py::extract<double>(p_art[t][1]));
	}
}

void PyInterface::init_adult_art_allocation(const double weight) {
	proj->dat.art_mort_weight(weight);
}

void PyInterface::init_adult_art_interruption(np::ndarray& art_exit_prop) {
	for (int t(0); t < proj->dat.num_years(); ++t) {
		proj->dat.art_exit_adult(t, DP::MALE,   py::extract<double>(art_exit_prop[t][0]));
		proj->dat.art_exit_adult(t, DP::FEMALE, py::extract<double>(art_exit_prop[t][0]));
	}
}

void PyInterface::init_adult_art_suppressed(np::ndarray& art_supp_pct) {
	const int num_ages = 4;
	int col_m, col_f;
	for (int t(0); t < proj->dat.num_years(); ++t) {
		for (int a(0); a < DP::N_AGE_ADULT; ++a) {
			col_m = std::min(a / 10, num_ages - 1);
			col_f = col_m + num_ages;
			proj->dat.art_suppressed_adult(t, DP::MALE,   a, py::extract<double>(art_supp_pct[t][col_m]));
			proj->dat.art_suppressed_adult(t, DP::FEMALE, a, py::extract<double>(art_supp_pct[t][col_f]));
		}
	}
}

// Test here. If this works well, move into DPUtil
void PyInterface::init_male_circumcision_uptake(np::ndarray& uptake) {
	using boost::math::interpolators::pchip;

	const size_t n(17); // number of 5-year age groups
	std::vector<double> x{0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85};
	std::vector<double> y(n + 1);
	double rate, prop;
	double dy_bgn, dy_end; // derivatives at left and right boundaries

	y[0] = 0.0;
	for (int t(0); t < proj->dat.num_years(); ++t) {
		// Calculate cumulative exposure to circumcision uptake at the
		// boundaries of five-year age groups
		for (int a(0); a < n; ++a) {
			prop = 0.01 * py::extract<double>(uptake[t][a]);
			rate = -5.0 * log(1.0 - prop);
			y[a+1] = y[a] + rate;
		}

		// Boundary calculations below are designed to match signal::pchip in
		// R, which uses public Fortran code by one of the authors of PCHIP
		// (doi:10.1137/0717021). R produces different results than boost::pchip
		// does by default because R approximates boundary derivatives using a
		// quadratic formula while boost::pchip takes uses a linear formula.
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

void PyInterface::init_clhiv_agein(np::ndarray& clhiv) {
	boost::multi_array_ref<double, 2>* ptr_clhiv = generate2d<double>(clhiv);
	DP::set_clhiv_agein(proj->dat, *ptr_clhiv);
	delete ptr_clhiv;
}

void PyInterface::init_effect_vmmc(const double effect) {
	proj->dat.effect_vmmc(effect);
}

void PyInterface::init_effect_condom(const double effect) {
	proj->dat.effect_condom(effect);
}

void PyInterface::project(const int year_final) {
	proj->project(year_final);
}

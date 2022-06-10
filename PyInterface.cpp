#include "PyInterface.h"
#include <fstream>
#include <iostream>
#include <boost/math/interpolators/pchip.hpp>

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

void PyInterface::init_mean_duration_pwid(const double years_female, const double years_male) {
	DP::set_mean_keypop_duration(proj->dat, DP::FEMALE, DP::POP_PWID, years_female);
	DP::set_mean_keypop_duration(proj->dat, DP::MALE,   DP::POP_PWID, years_male);
}

void PyInterface::init_mean_duration_fsw(const double years) {
	DP::set_mean_keypop_duration(proj->dat, DP::FEMALE, DP::POP_FSW, years);
}

void PyInterface::init_mean_duration_msm(const double years) {
	DP::set_mean_keypop_duration(proj->dat, DP::MALE, DP::POP_MSM, years);
}

void PyInterface::init_size_pwid(const double prop_female, const double prop_male) {
	proj->dat.keypop_size(DP::FEMALE, DP::POP_PWID, prop_female);
	proj->dat.keypop_size(DP::MALE,   DP::POP_PWID, prop_male);
}

void PyInterface::init_size_fsw(const double prop) {
	proj->dat.keypop_size(DP::FEMALE, DP::POP_FSW, prop);
}

void PyInterface::init_size_msm(const double prop) {
	proj->dat.keypop_size(DP::MALE, DP::POP_MSM, prop);
}

void PyInterface::init_size_trans(const double prop_female, const double prop_male) {
	proj->dat.keypop_size(DP::FEMALE, DP::POP_TRANS, prop_female);
	proj->dat.keypop_size(DP::MALE,   DP::POP_TRANS, prop_male);
}

void PyInterface::init_age_pwid(const double loc_female, const double shp_female, const double loc_male, const double shp_male) {
	DP::set_keypop_age(proj->dat, DP::FEMALE, DP::POP_PWID, loc_female, shp_female);
	DP::set_keypop_age(proj->dat, DP::MALE,   DP::POP_PWID, loc_male,   shp_male);
}

void PyInterface::init_age_fsw(const double loc, const double shp) {
	DP::set_keypop_age(proj->dat, DP::FEMALE, DP::POP_FSW, loc, shp);
}

void PyInterface::init_age_msm(const double loc, const double shp) {
	DP::set_keypop_age(proj->dat, DP::MALE, DP::POP_MSM, loc, shp);
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

void PyInterface::init_adult_art_curr(np::ndarray& art_num, np::ndarray& art_pct) {
	for (int t(0); t < proj->dat.num_years(); ++t) {
		proj->dat.art_num_adult(t, DP::MALE,   py::extract<double>(art_num[t][0]));
		proj->dat.art_num_adult(t, DP::FEMALE, py::extract<double>(art_num[t][1]));
		proj->dat.art_perc_adult(t, DP::MALE,   0.01 * py::extract<double>(art_pct[t][0]));
		proj->dat.art_perc_adult(t, DP::FEMALE, 0.01 * py::extract<double>(art_pct[t][1]));
	}
}

void PyInterface::init_adult_art_dropout(np::ndarray& art_drop_pct) {
	for (int t(0); t < proj->dat.num_years(); ++t) {
		proj->dat.art_drop_adult(t, DP::MALE,   -log(1.0 - py::extract<double>(art_drop_pct[t][0])));
		proj->dat.art_drop_adult(t, DP::FEMALE, -log(1.0 - py::extract<double>(art_drop_pct[t][0])));
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

#include "PyInterface.h"
#include <iostream>

// Helper function. Returns a pointer to a two-dimensional boost::multi_array_ref
// wrapper around the data stored in a numpy ndarray. NumType must be implicitly
// castable from the arr dtype, otherwise this will probably cause some astonishing
// errors. If you call this, you own the multi_array_ref pointer. Use it responsibly
// and clean up after you are done with it.
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

	// TODO:
	// Throw exception if
	// * sizeof(arr.dtype()) != sizeof(NumType) [multi_array_ref datatype incompatible with arr dtype]
	// ptr=nullptr case hit [mauled beyond recognition?]
	// arr.get_nd() > 2 [array dimension is larger than multiarray dimension]
	// other arr.get_flags() values that might cause issues?

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

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

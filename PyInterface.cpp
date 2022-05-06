#include "PyInterface.h"
#include <iostream>

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
	// TODO: throw exceptions if the dimension, size, or arrangement of pasfrs5y is
	// insufficient (https://docs.python.org/3/extending/extending.html) or
	// dat allocation fails
	const int nrows(pasfrs5y.shape(0));
	const int ncols(pasfrs5y.shape(1));
	double* data_ptr = reinterpret_cast<double*>(pasfrs5y.get_data());
	DP::year_age_ref_t* dat;

	if (pasfrs5y.get_flags() & np::ndarray::C_CONTIGUOUS) {
		dat = new DP::year_age_ref_t(data_ptr, boost::extents[nrows][ncols], boost::c_storage_order());
	} else if (pasfrs5y.get_flags() & np::ndarray::F_CONTIGUOUS) {
		dat = new DP::year_age_ref_t(data_ptr, boost::extents[nrows][ncols], boost::fortran_storage_order());
	} else {
		std::cerr << "Panic!" << std::endl;
		dat = nullptr;
	}

	proj->dat.init_pasfrs_from_5yr(*dat);

	delete dat;
}

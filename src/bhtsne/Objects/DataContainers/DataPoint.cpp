/*
 * DataPoint.cpp
 *
 *  Created on: Mar 22, 2017
 *      Author: nick
 */

#include "DataPoint.hpp"

namespace bhtsne {

DataPoint::DataPoint() {
	D_ = 1;
	ind_ = std::numeric_limits<uint32_t>::max();
	x_ = NULL;
}
DataPoint::DataPoint(int D, int ind, double* x) {
	D_ = D;
	ind_ = ind;
	x_ = (double*) malloc(D_ * sizeof(double));
	for (int d = 0; d < D_; ++d) {
		x_[d] = x[d];
	}
}
DataPoint::DataPoint(const DataPoint& other) {
	// this makes a deep copy -- should not free anything
	if (this != &other) {
		D_ = other.dimensionality();
		ind_ = other.index();
		x_ = (double*) malloc(D_ * sizeof(double));
		for (int d = 0; d < D_; ++d) {
			x_[d] = other.x(d);
		}
	}
}

DataPoint::~DataPoint() {
	if (x_ != NULL) {
		free(x_);
	}
}

DataPoint& DataPoint::operator=(const DataPoint& other) {
	// assignment should free old object
	if (this != &other) {
		if (x_ != NULL){
			free(x_);
		}
		D_ = other.dimensionality();
		ind_ = other.index();
		x_ = (double*) malloc(D_ * sizeof(double));
		for (int d = 0; d < D_; ++d) {
			x_[d] = other.x(d);
		}
	}
	return *this;
}

uint32_t DataPoint::index() const {
	return ind_;
}
uint32_t DataPoint::dimensionality() const {
	return D_;
}
double DataPoint::x(uint32_t d) const {
	return x_[d];
}


}  // namespace bhtsne



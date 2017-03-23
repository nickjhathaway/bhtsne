#pragma once
/*
 * DataPoint.hpp
 *
 *  Created on: Mar 22, 2017
 *      Author: nick
 */

#include "bhtsne/common.h"

namespace bhtsne {


class DataPoint {
	uint32_t ind_;

public:
	double* x_;
	uint32_t D_;

	DataPoint();
	DataPoint(int D, int ind, double* x);
	DataPoint(const DataPoint& other) ;

	~DataPoint();

	DataPoint& operator=(const DataPoint& other);

	uint32_t index() const;
	uint32_t dimensionality() const ;
	double x(uint32_t d) const;
};

inline double euclidean_distance(const DataPoint &t1, const DataPoint &t2) {
	double dd = .0;
	double* x1 = t1.x_;
	double* x2 = t2.x_;
	double diff;
	for (int d = 0; d < t1.D_; ++d) {
		diff = (x1[d] - x2[d]);
		dd += diff * diff;
	}
	return sqrt(dd);
}

}  // namespace bhtsne




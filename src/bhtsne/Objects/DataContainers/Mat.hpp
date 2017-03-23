#pragma once
/*
 * Mat.hpp
 *
 *  Created on: Mar 22, 2017
 *      Author: nick
 */

#include "bhtsne/common.h"


namespace bhtsne {

class Mat {
public:
	const size_t n_rows;
	const size_t n_cols;
	const size_t num_elements_;
	std::shared_ptr<double> data_;
	double* ptr_;

	Mat(const size_t n_rows, const size_t n_cols);
	Mat(const std::vector<std::vector<double>>& m);

	void write(std::ofstream & out) const;

};

}  // namespace bhtsne



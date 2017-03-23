/*
 * Mat.cpp
 *
 *  Created on: Mar 22, 2017
 *      Author: nick
 */



#include "Mat.hpp"

#include "bhtsne/Objects/DataContainers/MallocPtr.hpp"



namespace bhtsne {

Mat::Mat(const size_t n_rows, const size_t n_cols) :
		n_rows(n_rows), n_cols(n_cols), num_elements_(n_rows * n_cols) {
	data_ = MallocPtr<double>::NumElements(num_elements_, true);
	ptr_ = data_.get();
}

Mat::Mat(const std::vector<std::vector<double>>& m) :
		n_rows(m.size()), n_cols(m.front().size()), num_elements_(n_rows * n_cols) {
	data_ = MallocPtr<double>::NumElements(num_elements_, false);
	ptr_ = data_.get();

	for (size_t i = 0; i < n_rows; ++i) {
		const auto& p = m[i];
		std::copy(p.begin(), p.end(), &ptr_[i * n_cols]);
	}
}

void Mat::write(std::ofstream & out) const{
	for(uint32_t i = 0; i < num_elements_; ++i){
		uint32_t colPos = i % n_cols;
		//uint32_t rowPos = i / output.n_cols;
		if(0 != colPos){
			out << "\t";
		}
		out << ptr_[i];
		if(colPos + 1 == n_cols){
			out << "\n";
		}
	}
	//make sure everything is written to file and not still in buffer
	//just a precaution against someone trying to read from out before it closes
	out.flush();
}

}  // namespace bhtsne

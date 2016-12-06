#include <cppitertools/range.hpp>
#include <sstream>
#include <ostream>
#include <istream>
#include <iostream>
#include <string>
#include <armadillo>
#include <pca.h>

#include "tsne.h"

Mat loadData(const TSNEArgs& params){
  std::ifstream inFile("data.tsv");

  if(!inFile){
    std::cerr << "Error in reading in data.tsv" << std::endl;
    exit(1);
  }

  std::string line = "";
  uint32_t lineCount = 0;
  std::vector<std::vector<double>> inputVecVec;
  while(std::getline(inFile, line)){
    ++lineCount;
    std::cout << "\r"<< lineCount;
    std::cout.flush();
    std::vector<double> tempVec;
    std::stringstream ss;
    ss << line;
    double out = 0;
    while(ss >> out){
      tempVec.emplace_back(out);
    }
    if(!inputVecVec.empty()){
    	if(inputVecVec.front().size() != tempVec.size()){
    		std::stringstream ss;
    		 ss << __PRETTY_FUNCTION__ << ", error input vector is not the size as the rest of the matrix" << "\n";
    		 ss << "Matrix size: " << inputVecVec.front().size() << "\n";
    		 ss << "Adding size: " << tempVec.size() << "\n";
    		 throw std::runtime_error{ss.str()};
    	}
    }
    inputVecVec.emplace_back(tempVec);
  }
  std::cout << std::endl;

  uint32_t colNum = inputVecVec.front().size();
  uint32_t rowNum = inputVecVec.size();
  if(params.menaCenterCols_){
    for(const auto & colPos : iter::range(colNum)){
    	double sum = 0;
    	for(const auto & rowPos : iter::range(inputVecVec.size())){
    		sum += inputVecVec[rowPos][colPos];
    	}
    	double mean = sum/rowNum;
    	for(const auto & rowPos : iter::range(inputVecVec.size())){
    		inputVecVec[rowPos][colPos] = inputVecVec[rowPos][colPos] - mean;
    	}
    }
  }


	std::vector<std::vector<double>> princomp;
  if(params.doPca_){
  	stats::pca pca(inputVecVec.front().size());
  	for (const auto & row : inputVecVec) {
  		pca.add_record(row);
  	}
  	pca.set_solver("standard");
  	pca.solve();
  	pca.set_num_retained(params.initial_dim_);


  	for(const auto & row : inputVecVec){
  		princomp.emplace_back(pca.to_principal_space(row));
  	}
  	/*
  	std::ofstream outPca("temp_out_pca.tsv");
    for(const auto  rowPos : iter::range(princomp.size())){
      for(const auto colPos : iter::range(princomp[rowPos].size())){
        if(colPos!= 0){
        	outPca << "\t";
        }
        outPca << princomp[rowPos][colPos];
      }
      outPca << std::endl;
    }*/
  }else{
  	princomp = inputVecVec;
  }
  Mat input(princomp);
  return input;
}

int main() {
	TSNEArgs params;
	params.perplexity_ = 50;
	const auto input = loadData(params);

	TSNE trunner;
	auto output = trunner.run(input, params);

	std::vector<std::vector<double>> trueOutput(output.n_rows, std::vector<double>(output.n_cols, 0));
	std::ofstream outFile("temp_out.txt");
	//armadillo is column wise so the output has to go down the columns
	uint32_t count = 0;
	/*
	for (const auto colPos : iter::range(output.n_cols)) {
		for (const auto rowPos : iter::range(output.n_rows)) {
			trueOutput[count / output.n_cols][count % output.n_cols] = output(rowPos,
					colPos);
			++count;
		}
	}

	for (const auto rowPos : iter::range(trueOutput.size())) {
		for (const auto colPos : iter::range(trueOutput[rowPos].size())) {
			outFile << trueOutput[rowPos][colPos];
			if (colPos + 1 != trueOutput[rowPos].size()) {
				outFile << "\t";
			}
		}
		outFile << std::endl;
	}
	*/
}



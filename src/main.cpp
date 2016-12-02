
#include "tsne.h"
#include <cppitertools/range.hpp>
#include <sstream>
#include <ostream>
#include <istream>
#include <iostream>
#include <string>
#include <pca.h>

arma::mat loadData(const TSNEArgs& params){
	bool doPca = true;
	bool menaCenterCols = true;

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
  if(menaCenterCols){
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
  if(doPca){
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
  	std::ofstream outPca("temp_out_pca.tsv");
    for(const auto  rowPos : iter::range(princomp.size())){
      for(const auto colPos : iter::range(princomp[rowPos].size())){
        if(colPos!= 0){
        	outPca << "\t";
        }
        outPca << princomp[rowPos][colPos];
      }
      outPca << std::endl;
    }
  }else{
  	princomp = inputVecVec;
  }
  arma::mat input(princomp.size(), princomp.front().size());

  for(const auto  rowPos : iter::range(princomp.size())){
    for(const auto colPos : iter::range(princomp[rowPos].size())){
      input(rowPos, colPos) = princomp[rowPos][colPos];
    }
  }
  return input;
}

int main() {
	TSNEArgs params;
	//params.max_iter_ = 400;
	params.perplexity_ = 50;
	const auto input = loadData(params);

	TSNE trunner;
	auto output = trunner.run(input, params);

	std::ofstream outFile("temp_outOut.txt");
	for (const auto rowPos : iter::range(output.n_rows)) {
		for (const auto colPos : iter::range(output.n_cols)) {
			outFile << output(rowPos, colPos);
			if (colPos + 1 != output.n_cols) {
				outFile << "\t";
			}
		}
		outFile << std::endl;
	}

}

/*
using namespace std;

int test() {

	const int num_variables = 10;
	const int num_records = 300;

	stats::pca pca(num_variables);
	pca.set_do_bootstrap(true, 100);

	cout<<"Adding random data records ..."<<endl;
	srand(1);
	for (int i=0; i<num_records; ++i) {
		vector<double> record(num_variables);
		for (auto& value : record) {
			value = rand()%20 - 10;
		}
		pca.add_record(record);
	}

	cout<<"Solving ..."<<endl;
	pca.solve();

	cout<<"Energy = "<<pca.get_energy()<<" ("<<
      		stats::utils::get_sigma(pca.get_energy_boot())<<")"<<endl;

	const auto eigenvalues = pca.get_eigenvalues();
	cout<<"First three eigenvalues = "<<eigenvalues[0]<<", "
									  <<eigenvalues[1]<<", "
									  <<eigenvalues[2]<<endl;

	cout<<"Orthogonal Check = "<<pca.check_eigenvectors_orthogonal()<<endl;
	cout<<"Projection Check = "<<pca.check_projection_accurate()<<endl;

	pca.save("pca_results");

	return 0;
}*/

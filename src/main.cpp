
#include "tsne.h"
#include <cppitertools/range.hpp>
#include <sstream>
#include <ostream>
#include <istream>
#include <iostream>
#include <string>



int main(){
	
	TSNEArgs params;

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
		inputVecVec.emplace_back(tempVec);
	}
	std::cout << std::endl;
	std::cout << "inputVecVec.front().size(): " << inputVecVec.front().size() << std::endl;
	arma::mat input(inputVecVec.size(), inputVecVec.front().size());
	std::cout << input.n_rows << std::endl;
	std::cout << input.n_cols << std::endl;
	for(const auto  rowPos : iter::range(inputVecVec.size())){
		for(const auto colPos : iter::range(inputVecVec[rowPos].size())){
			input(rowPos, colPos) = inputVecVec[rowPos][colPos];
		}
	}
	TSNE trunner;

//params.max_iter_ = 400;
	params.perplexity_ = 150;
	auto output = trunner.run(input, params);
	std::ofstream outFile("temp_outOut.txt");
	for(const auto  rowPos : iter::range(output.n_rows)){
		for(const auto colPos : iter::range(output.n_cols)){
			outFile << output(rowPos, colPos);
			if(colPos + 1 != output.n_cols){
				outFile<< "\t";
			}
		}
		outFile <<std::endl;
	}

}

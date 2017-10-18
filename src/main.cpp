#include <cppitertools/range.hpp>
#include <sstream>
#include <ostream>
#include <istream>
#include <iostream>
#include <string>
#include <armadillo>
#include <pca.h>

#include "bhtsne.h"
#include <bibcpp/progutils.h>


namespace bhtsne {


void meanCenterColsMatrix(std::vector<std::vector<double>> & mat) {
	uint32_t colNum = mat.front().size();
	uint32_t rowNum = mat.size();
	for (const auto & colPos : iter::range(colNum)) {
		double sum = 0;
		for (const auto & rowPos : iter::range(mat.size())) {
			sum += mat[rowPos][colPos];
		}
		double mean = sum / rowNum;
		for (const auto & rowPos : iter::range(mat.size())) {
			mat[rowPos][colPos] = mat[rowPos][colPos] - mean;
		}
	}
}


struct pcaArgs{
	enum class SolverType{
		STANDARD,
		DC
	};
	SolverType type_ = SolverType::STANDARD;
	uint32_t numRetained_ = std::numeric_limits<uint32_t>::max(); /**< Number of components to retain */
};


std::vector<std::vector<double>> runPca(
		const std::vector<std::vector<double>> & input, pcaArgs args) {
	std::vector<std::vector<double>> ret;
	if(0 == input.size()){
		//maybe throw?
		return ret;
	}
	stats::pca pca(input.front().size());
	for (const auto & row : input) {
		pca.add_record(row);
	}
	switch (args.type_) {
	case pcaArgs::SolverType::DC:
		pca.set_solver("dc");
		break;
	case pcaArgs::SolverType::STANDARD:
		pca.set_solver("standard");
		break;
	default:
		break;
	}

	pca.solve();
	if (std::numeric_limits<uint32_t>::max() != args.numRetained_
			&& args.numRetained_ < input.front().size()) {
		pca.set_num_retained(args.numRetained_);
	}
	for (const auto & row : input) {
		ret.emplace_back(pca.to_principal_space(row));
	}
	return ret;
}

std::vector<std::vector<double>> readInData(const bib::bfs::path & inputFnp, bool verbose){
	if(!bib::bfs::exists(inputFnp)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, input file: " << inputFnp << " doesn't exist " << "\n";
		throw std::runtime_error{ss.str()};
	}
  std::string line = "";
  uint32_t lineCount = 0;
  std::vector<std::vector<double>> ret;
  std::ifstream inFile(inputFnp.string());
  if(!inFile){
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, in opening file: " << inputFnp << "\n";
		throw std::runtime_error{ss.str()};
  }
  while(bib::files::crossPlatGetline(inFile, line)){
    ++lineCount;
    if(verbose){
      std::cout << "\r"<< lineCount;
      std::cout.flush();
    }
    std::vector<double> tempVec;
    std::stringstream ss;
    ss << line;
    double out = 0;
    while(ss >> out){
      tempVec.emplace_back(out);
    }
    if(!ret.empty()){
    	if(ret.front().size() != tempVec.size()){
    		std::stringstream ss;
    		 ss << __PRETTY_FUNCTION__ << ", error input vector is not the size as the rest of the matrix" << "\n";
    		 ss << "Matrix size: " << ret.front().size() << "\n";
    		 ss << "Adding size: " << tempVec.size() << "\n";
    		 throw std::runtime_error{ss.str()};
    	}
    }
    ret.emplace_back(tempVec);
  }
  if(verbose){
  	std::cout << std::endl;
  }
  return ret;
}

Mat preProcessData(std::vector<std::vector<double>> & inputVecVec, const TSNEArgs& params){
  std::ifstream inFile("data.tsv");

  if(params.meanCenterCols_){
  	meanCenterColsMatrix(inputVecVec);
  }
	std::vector<std::vector<double>> princomp;
  if(params.doPca_){
  	pcaArgs pArgs;
  	pArgs.numRetained_ = params.initial_dim_;
  	princomp = runPca(inputVecVec, pArgs);
  }else{
  	princomp = inputVecVec;
  }
  Mat input(princomp);
  return input;
}

Mat loadData(const bib::bfs::path & inputFnp, const TSNEArgs& params){
  //std::ifstream inFile("data.tsv");
  auto inputVecVec = readInData(inputFnp, params.verbose_);
  return preProcessData(inputVecVec, params);
}



}  // namespace bhtsne

int main(int argc, char* argv[]) {
	bhtsne::TSNEArgs params;
	params.perplexity_ = 50;
	params.no_dims_ = 3;
	bib::bfs::path inputFile = "";
	bib::bfs::path outputFile = "out.tab.txt";
	bool overWrite = false;
	bool noPca = false;
	bib::progutils::ProgramSetUp setUp(argc, argv);
	setUp.setOption(params.perplexity_, "--perplexity", "TSNE perplexity");
	setUp.setOption(params.no_dims_, "--noDims", "TSNE output dimensions");
	setUp.setOption(params.initial_dim_, "--initialDims", "Initial input dimensions");
	setUp.setOption(params.verbose_, "--verbose", "Run verbose");

	setUp.setOption(outputFile, "--out", "Output file name");
	setUp.setOption(noPca, "--noPca", "Don't do pca");
	params.doPca_= !noPca;
	setUp.setOption(inputFile, "--in", "Input file name, should be a whitespace delimited matrix", true);
	setUp.setOption(overWrite, "--overWrite", "Overwrite output file if it exists");

	setUp.finishSetUp(std::cout);

	//check if ouput file exists already
	if(bib::bfs::exists(outputFile) && !overWrite){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, file " << outputFile << " already exists, use --overWrite to over write it" << "\n";
		throw std::runtime_error{ss.str()};
	}

	//read in and process data
	const auto input = loadData(inputFile, params);
	//run tsne
	bhtsne::TSNE trunner;
	auto output = trunner.run(input, params);
	//write output
	std::ofstream outFile(outputFile.string());
	output.write(outFile);
	outFile.close();
}



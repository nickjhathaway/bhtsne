
#include "tsne.h"
#include <cppitertools/range.hpp>
#include <sstream>
#include <ostream>
#include <istream>
#include <iostream>
#include <string>
#include <pca.h>

arma::mat loadData(){
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

  return input;
}

int main(){
  const auto input = loadData();
  
  TSNEArgs params;
  //params.max_iter_ = 400;
  params.perplexity_ = 50;
	
  TSNE trunner;
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
}

#pragma once
/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

#define ARMA_64BIT_WORD
#include <armadillo>

struct TSNEArgs{

  double theta_ = 0.5;        // gradient accuracy
  double perplexity_ = 30;     // perplexity
  uint32_t no_dims_ = 2;         // output dimensionality
  bool skip_random_init_ = false;
  uint32_t max_iter_ = 1000;    // maximum number of iterations
  uint32_t stop_lying_iter_=250;
  uint32_t mom_switch_iter_=250;
  uint32_t rand_seed_ = -1;

};

static inline double sign(double x) { return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0)); }

class TSNE {
public:
  arma::mat run(const arma::mat& X, const TSNEArgs& params);
  
  void run(double* X, uint32_t N, uint32_t D, double* Y, uint32_t no_dims, double perplexity,
	     double theta, int rand_seed,
             bool skip_random_init, uint32_t max_iter=1000, uint32_t stop_lying_iter=250,
	     uint32_t mom_switch_iter=250);
    
    bool load_data(double** data, int* n, int* d, int* no_dims, double* theta,
		   double* perplexity, int* rand_seed, int* max_iter);
    void save_data(double* data, int* landmarks, double* costs, uint32_t n, uint32_t d);
    void symmetrizeMatrix(uint32_t** row_P, uint32_t** col_P,
			  double** val_P, uint32_t N); // should be static!

private:
    void computeGradient(double* P, uint32_t* inp_row_P, uint32_t* inp_col_P,
			 double* inp_val_P, double* Y, uint32_t N, uint32_t D, double* dC, double theta);
    void computeExactGradient(double* P, double* Y, uint32_t N, uint32_t D, double* dC);
    double evaluateError(double* P, double* Y, uint32_t N, uint32_t D);
    double evaluateError(uint32_t* row_P, uint32_t* col_P, double* val_P, double* Y,
			 uint32_t N, uint32_t D, double theta);
    void zeroMean(double* X, uint32_t N, uint32_t D);
    void computeGaussianPerplexity(double* X, uint32_t N, uint32_t D, double* P, double perplexity);
    void computeGaussianPerplexity(double* X, uint32_t N, uint32_t D, uint32_t** _row_P, uint32_t** _col_P, double** _val_P, double perplexity, uint32_t K);
    void computeSquaredEuclideanDistance(double* X, uint32_t N, uint32_t D, double* DD);
    double randn();
};

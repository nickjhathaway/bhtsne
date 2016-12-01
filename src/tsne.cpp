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



#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>
#include "vptree.h"
#include "sptree.h"
#include "tsne.h"

#include <cppitertools/range.hpp>

using namespace std;

// Perform t-SNE
arma::mat TSNE::run(const arma::mat& X, const TSNEArgs& params) {
  arma::mat coeff;
  arma::mat score;
  arma::vec latent;
  arma::vec tsquared;
  arma::princomp(coeff, score, latent, tsquared, X);

  arma::mat x = X * coeff;
  arma::mat ret(x.n_rows, params.no_dims_);

  run(x.memptr(), x.n_rows, x.n_cols, ret.memptr(),
      params.no_dims_, params.perplexity_, params.theta_, params.rand_seed_,
      params.skip_random_init_, params.max_iter_, params.stop_lying_iter_,
      params.mom_switch_iter_);
  
  return ret;
}

// Perform t-SNE
void TSNE::run(double* X, uint32_t N, uint32_t D, double* Y,
	       uint32_t no_dims, double perplexity, double theta, int rand_seed,
               bool skip_random_init, uint32_t max_iter, uint32_t stop_lying_iter,
	       uint32_t mom_switch_iter) {

  // Set random seed
  if (skip_random_init != true) {
    if(rand_seed >= 0) {
      printf("Using random seed: %d\n", rand_seed);
      srand((uint32_t) rand_seed);
    } else {
      printf("Using current time as random seed...\n");
      srand(time(NULL));
    }
  }

  // Determine whether we are using an exact algorithm
  std::cout << "N" << N << std::endl;
  if(N - 1 < 3 * perplexity) {
    printf("Perplexity too large for the number of data points!\n");
    exit(1);
  }
  printf("Using no_dims = %d, perplexity = %f, and theta = %f\n", no_dims, perplexity, theta);
  bool exact = (theta == .0) ? true : false;

  // Set learning parameters
  float total_time = .0;
  clock_t start, end;
  double momentum = .5, final_momentum = .8;
  double eta = 200.0;

  // Allocate some memory
  double* dY    = (double*) malloc(N * no_dims * sizeof(double));
  double* uY    = (double*) malloc(N * no_dims * sizeof(double));
  double* gains = (double*) malloc(N * no_dims * sizeof(double));
  if(dY == nullptr || uY == nullptr || gains == nullptr) {
    printf("Memory allocation failed!\n");
    exit(1);
  }

  for(uint32_t i = 0; i < N * no_dims; i++){
    uY[i] =  .0;
  }
  for(uint32_t i = 0; i < N * no_dims; i++) {
    gains[i] = 1.0;
  }

  // Normalize input data (to prevent numerical problems)
  printf("Computing input similarities...\n");
  start = clock();
  zeroMean(X, N, D);
  double max_X = .0;

  for(uint32_t i = 0; i < N * D; i++) {
    if(fabs(X[i]) > max_X) {
      max_X = fabs(X[i]);
    }
  }

  for(uint32_t i = 0; i < N * D; i++) {
    X[i] /= max_X;
  }

  // Compute input similarities for exact t-SNE
  double* P = nullptr;
  uint32_t* row_P = nullptr;
  uint32_t* col_P = nullptr;
  double* val_P = nullptr;

  if(exact) {

    // Compute similarities
    printf("Exact?");
    P = (double*) malloc(N * N * sizeof(double));
    if(!P) {
      printf("Memory allocation failed!\n");
      exit(1);
    }
    computeGaussianPerplexity(X, N, D, P, perplexity);

    // Symmetrize input similarities
    printf("Symmetrizing...\n");
    uint32_t nN = 0;
    for (uint32_t n = 0; n < N; n++) {
      uint32_t mN = (n + 1) * N;
      for (uint32_t m = n + 1; m < N; m++) {
	P[nN + m] += P[mN + n];
	P[mN + n] = P[nN + m];
	mN += N;
      }
      nN += N;
    }
    double sum_P = .0;
    for(uint32_t i = 0; i < N * N; i++) {
      sum_P += P[i];
    }
    for(uint32_t i = 0; i < N * N; i++) {
      P[i] /= sum_P;
    }
  } else {
    // Compute input similarities for approximate t-SNE
    
    // Compute asymmetric pairwise input similarities
    computeGaussianPerplexity(X, N, D, &row_P, &col_P, &val_P, perplexity,
			      uint32_t{3 * perplexity});

    // Symmetrize input similarities
    symmetrizeMatrix(&row_P, &col_P, &val_P, N);
    double sum_P = .0;
    for(uint32_t i = 0; i < row_P[N]; i++) {
      sum_P += val_P[i];
    }
    for(uint32_t i = 0; i < row_P[N]; i++) {
      val_P[i] /= sum_P;
    }
  }
  end = clock();

  // Lie about the P-values
  if(exact) {
    for(uint32_t i = 0; i < N * N; i++){
      P[i] *= 12.0;
    }
  } else {
    for(uint32_t i = 0; i < row_P[N]; i++) {
      val_P[i] *= 12.0;
    }
  }

  // Initialize solution (randomly)
  if (skip_random_init != true) {
    for(uint32_t i = 0; i < N * no_dims; i++) {
      Y[i] = randn() * .0001;
    }
  }

  // Perform main training loop
  if(exact) {
    printf("Input similarities computed in %4.2f seconds!\nLearning embedding...\n",
	   (float) (end - start) / CLOCKS_PER_SEC);
  } else {
    printf("Input similarities computed in %4.2f seconds (sparsity = %f)!\nLearning embedding...\n",
	   (float) (end - start) / CLOCKS_PER_SEC, (double) row_P[N] / ((double) N * (double) N));
  }
  start = clock();

  for(uint32_t iter = 0; iter < max_iter; ++iter) {

    // Compute (approximate) gradient
    if(exact) {
      computeExactGradient(P, Y, N, no_dims, dY);
    } else {
      computeGradient(P, row_P, col_P, val_P, Y, N, no_dims, dY, theta);
    }

    // Update gains
    for(uint32_t i = 0; i < N * no_dims; i++) {
      gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
    }
    for(uint32_t i = 0; i < N * no_dims; i++) {
      if(gains[i] < .01) gains[i] = .01;
    }

    // Perform gradient update (with momentum and gains)
    for(uint32_t i = 0; i < N * no_dims; i++) {
      uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
    }
    for(uint32_t i = 0; i < N * no_dims; i++) {
      Y[i] = Y[i] + uY[i];
    }

    // Make solution zero-mean
    zeroMean(Y, N, no_dims);

    // Stop lying about the P-values after a while, and switch momentum
    if (iter == stop_lying_iter) {
      if (exact) {
	for (uint32_t i = 0; i < N * N; i++) {
	  P[i] /= 12.0;
	}
      } else {
	for (uint32_t i = 0; i < row_P[N]; i++) {
	  val_P[i] /= 12.0;
	}
      }
    }
    if (iter == mom_switch_iter) {
      momentum = final_momentum;
    }

    // Print out progress
    if (iter > 0 && (iter % 50 == 0 || iter == max_iter - 1)) {
      end = clock();
      double C = .0;
      if (exact) {
	C = evaluateError(P, Y, N, no_dims);
      } else {
	C = evaluateError(row_P, col_P, val_P, Y, N, no_dims, theta); // doing approximate computation here!
      }

      if (iter == 0) {
	printf("Iteration %d: error is %f\n", iter + 1, C);
      } else {
	total_time += (float) (end - start) / CLOCKS_PER_SEC;
	printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n",
	       iter, C, (float) (end - start) / CLOCKS_PER_SEC);
      }
      start = clock();
    }
  }
  end = clock(); total_time += (float) (end - start) / CLOCKS_PER_SEC;

  // Clean up memory
  free(dY);
  free(uY);
  free(gains);
  if(exact) {
    free(P);
  } else {
    free(row_P); row_P = nullptr;
    free(col_P); col_P = nullptr;
    free(val_P); val_P = nullptr;
  }
  printf("Fitting performed in %4.2f seconds.\n", total_time);
}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
void TSNE::computeGradient(double* P, uint32_t* inp_row_P, uint32_t* inp_col_P, double* inp_val_P, double* Y, uint32_t N, uint32_t D, double* dC, double theta)
{

  // Construct space-partitioning tree on current map
  SPTree tree(D, Y, N);

  // Compute all terms required for t-SNE gradient
  double sum_Q = .0;
  double* pos_f = (double*) calloc(N * D, sizeof(double));
  double* neg_f = (double*) calloc(N * D, sizeof(double));
  if(pos_f == nullptr || neg_f == nullptr) { printf("Memory allocation failed!\n"); exit(1); }
  tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);
  
  for(uint32_t n = 0; n < N; n++) {
    tree.computeNonEdgeForces(n, theta, neg_f + n * D, &sum_Q);
  }

  // Compute final t-SNE gradient
  for(uint32_t i = 0; i < N * D; i++) {
    dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
  }
  free(pos_f);
  free(neg_f);
}

// Compute gradient of the t-SNE cost function (exact)
void TSNE::computeExactGradient(double* P, double* Y, uint32_t N, uint32_t D, double* dC) {

  // Make sure the current gradient contains zeros
  for(uint32_t i = 0; i < N * D; i++) {
    dC[i] = 0.0;
  }

  // Compute the squared Euclidean distance matrix
  double* DD = (double*) malloc(N * N * sizeof(double));
  if(DD == nullptr) { printf("Memory allocation failed!\n"); exit(1); }
  computeSquaredEuclideanDistance(Y, N, D, DD);

  // Compute Q-matrix and normalization sum
  double* Q    = (double*) malloc(N * N * sizeof(double));
  if(Q == nullptr) { printf("Memory allocation failed!\n"); exit(1); }
  double sum_Q = .0;
  uint32_t nN = 0;
  for(uint32_t n = 0; n < N; n++) {
    for(uint32_t m = 0; m < N; m++) {
      if(n != m) {
	Q[nN + m] = 1 / (1 + DD[nN + m]);
	sum_Q += Q[nN + m];
      }
    }
    nN += N;
  }

  // Perform the computation of the gradient
  nN = 0;
  uint32_t nD = 0;
  for(uint32_t n = 0; n < N; n++) {
    uint32_t mD = 0;
    for(uint32_t m = 0; m < N; m++) {
      if(n != m) {
	double mult = (P[nN + m] - (Q[nN + m] / sum_Q)) * Q[nN + m];
	for(uint32_t d = 0; d < D; d++) {
	  dC[nD + d] += (Y[nD + d] - Y[mD + d]) * mult;
	}
      }
      mD += D;
    }
    nN += N;
    nD += D;
  }

  // Free memory
  free(DD); DD = nullptr;
  free(Q);  Q  = nullptr;
}


// Evaluate t-SNE cost function (exactly)
double TSNE::evaluateError(double* P, double* Y, uint32_t N, uint32_t D) {

  // Compute the squared Euclidean distance matrix
  double* DD = (double*) malloc(N * N * sizeof(double));
  double* Q = (double*) malloc(N * N * sizeof(double));
  if(DD == nullptr || Q == nullptr) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  computeSquaredEuclideanDistance(Y, N, D, DD);

  // Compute Q-matrix and normalization sum
  uint32_t nN = 0;
  double sum_Q = DBL_MIN;
  for(uint32_t n = 0; n < N; n++) {
    for(uint32_t m = 0; m < N; m++) {
      if(n != m) {
	Q[nN + m] = 1 / (1 + DD[nN + m]);
	sum_Q += Q[nN + m];
      }
      else Q[nN + m] = DBL_MIN;
    }
    nN += N;
  }
  for(uint32_t i = 0; i < N * N; i++) {
    Q[i] /= sum_Q;
  }

  // Sum t-SNE error
  double C = .0;
  for(uint32_t n = 0; n < N * N; n++) {
    C += P[n] * log((P[n] + FLT_MIN) / (Q[n] + FLT_MIN));
  }

  // Clean up memory
  free(DD);
  free(Q);
  return C;
}

// Evaluate t-SNE cost function (approximately)
double TSNE::evaluateError(uint32_t* row_P, uint32_t* col_P, double* val_P, double* Y, uint32_t N, uint32_t D, double theta){
  // Get estimate of normalization term
  SPTree tree(D, Y, N);
  double* buff = (double*) calloc(D, sizeof(double));
  double sum_Q = .0;
  for(uint32_t n = 0; n < N; n++) {
    tree.computeNonEdgeForces(n, theta, buff, &sum_Q);
  }

  // Loop over all edges to compute t-SNE error
  uint32_t ind1, ind2;
  double C = .0, Q;
  for(uint32_t n = 0; n < N; n++) {
    ind1 = n * D;
    for(uint32_t i = row_P[n]; i < row_P[n + 1]; i++) {
      Q = .0;
      ind2 = col_P[i] * D;
      for(uint32_t d = 0; d < D; d++) buff[d]  = Y[ind1 + d];
      for(uint32_t d = 0; d < D; d++) buff[d] -= Y[ind2 + d];
      for(uint32_t d = 0; d < D; d++) Q += buff[d] * buff[d];
      Q = (1.0 / (1.0 + Q)) / sum_Q;
      C += val_P[i] * log((val_P[i] + FLT_MIN) / (Q + FLT_MIN));
    }
  }

  // Clean up memory
  free(buff);
  return C;
}


// Compute input similarities with a fixed perplexity
void TSNE::computeGaussianPerplexity(double* X, uint32_t N, uint32_t D, double* P, double perplexity) {

  // Compute the squared Euclidean distance matrix
  double* DD = (double*) malloc(N * N * sizeof(double));
  if(DD == nullptr) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  computeSquaredEuclideanDistance(X, N, D, DD);

  // Compute the Gaussian kernel row by row
  uint32_t nN = 0;
  for(uint32_t n = 0; n < N; n++) {

    // Initialize some variables
    bool found = false;
    double beta = 1.0;
    double min_beta = -DBL_MAX;
    double max_beta =  DBL_MAX;
    double tol = 1e-5;
    double sum_P;

    // Iterate until we found a good perplexity
    uint32_t iter = 0;
    while(!found && iter < 200) {

      // Compute Gaussian kernel row
      for(uint32_t m = 0; m < N; m++) P[nN + m] = exp(-beta * DD[nN + m]);
      P[nN + n] = DBL_MIN;

      // Compute entropy of current row
      sum_P = DBL_MIN;
      for(uint32_t m = 0; m < N; m++) sum_P += P[nN + m];
      double H = 0.0;
      for(uint32_t m = 0; m < N; m++) H += beta * (DD[nN + m] * P[nN + m]);
      H = (H / sum_P) + log(sum_P);

      // Evaluate whether the entropy is within the tolerance level
      double Hdiff = H - log(perplexity);
      if(Hdiff < tol && -Hdiff < tol) {
	found = true;
      }
      else {
	if(Hdiff > 0) {
	  min_beta = beta;
	  if(max_beta == DBL_MAX || max_beta == -DBL_MAX) {
	    beta *= 2.0;
	  } else {
	    beta = (beta + max_beta) / 2.0;
	  }
	}
	else {
	  max_beta = beta;
	  if(min_beta == -DBL_MAX || min_beta == DBL_MAX) {
	    beta /= 2.0;
	  } else {
	    beta = (beta + min_beta) / 2.0;
	  }
	}
      }

      // Update iteration counter
      iter++;
    }

    // Row normalize P
    for(uint32_t m = 0; m < N; m++) {
      P[nN + m] /= sum_P;
    }
    nN += N;
  }

  // Clean up memory
  free(DD); DD = nullptr;
}


// Compute input similarities with a fixed perplexity using ball trees (this function allocates memory another function should free)
void TSNE::computeGaussianPerplexity(double* X, uint32_t N, uint32_t D, uint32_t** _row_P,
				     uint32_t** _col_P, double** _val_P, double perplexity,
				     uint32_t K) {

  if(perplexity > K) {
    printf("Perplexity should be lower than K!\n");
  }

  // Allocate the memory we need
  *_row_P = (uint32_t*)    malloc((N + 1) * sizeof(uint32_t));
  *_col_P = (uint32_t*)    calloc(N * K, sizeof(uint32_t));
  *_val_P = (double*) calloc(N * K, sizeof(double));
  if(*_row_P == nullptr || *_col_P == nullptr || *_val_P == nullptr) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  uint32_t* row_P = *_row_P;
  uint32_t* col_P = *_col_P;
  double* val_P = *_val_P;
  double* cur_P = (double*) malloc((N - 1) * sizeof(double));
  if(cur_P == nullptr) {
    printf("Memory allocation failed!\n"); exit(1);
  }
  row_P[0] = 0;
  for(uint32_t n = 0; n < N; n++) {
    row_P[n + 1] = row_P[n] + (uint32_t) K;
  }

  // Build ball tree on data set
  VpTree<DataPoint, euclidean_distance> tree;
  vector<DataPoint> obj_X(N, DataPoint(D, -1, X));
  for(uint32_t n = 0; n < N; n++) {
    obj_X[n] = DataPoint(D, n, X + n * D);
  }
  tree.create(obj_X);

  // Loop over all points to find nearest neighbors
  printf("Building tree...\n");
  vector<DataPoint> indices;
  vector<double> distances;
  for(uint32_t n = 0; n < N; n++) {

    if(n % 10000 == 0) printf(" - point %d of %d\n", n, N);

    // Find nearest neighbors
    indices.clear();
    distances.clear();
    tree.search(obj_X[n], K + 1, &indices, &distances);

    // Initialize some variables for binary search
    bool found = false;
    double beta = 1.0;
    double min_beta = -DBL_MAX;
    double max_beta =  DBL_MAX;
    double tol = 1e-5;

    // Iterate until we found a good perplexity
    uint32_t iter = 0; double sum_P;
    while(!found && iter < 200) {

      // Compute Gaussian kernel row
      for(uint32_t m = 0; m < K; m++) cur_P[m] = exp(-beta * distances[m + 1] * distances[m + 1]);

      // Compute entropy of current row
      sum_P = DBL_MIN;
      for(uint32_t m = 0; m < K; m++) sum_P += cur_P[m];
      double H = .0;
      for(uint32_t m = 0; m < K; m++) H += beta * (distances[m + 1] * distances[m + 1] * cur_P[m]);
      H = (H / sum_P) + log(sum_P);

      // Evaluate whether the entropy is within the tolerance level
      double Hdiff = H - log(perplexity);
      if(Hdiff < tol && -Hdiff < tol) {
	found = true;
      }
      else {
	if(Hdiff > 0) {
	  min_beta = beta;
	  if(max_beta == DBL_MAX || max_beta == -DBL_MAX)
	    beta *= 2.0;
	  else
	    beta = (beta + max_beta) / 2.0;
	}
	else {
	  max_beta = beta;
	  if(min_beta == -DBL_MAX || min_beta == DBL_MAX)
	    beta /= 2.0;
	  else
	    beta = (beta + min_beta) / 2.0;
	}
      }

      // Update iteration counter
      iter++;
    }

    // Row-normalize current row of P and store in matrix
    for(uint32_t m = 0; m < K; m++) {
      cur_P[m] /= sum_P;
    }
    for(uint32_t m = 0; m < K; m++) {
      col_P[row_P[n] + m] = uint32_t{indices[m + 1].index()};
      val_P[row_P[n] + m] = cur_P[m];
    }
  }

  // Clean up memory
  obj_X.clear();
  free(cur_P);
}


// Symmetrizes a sparse matrix
void TSNE::symmetrizeMatrix(uint32_t** _row_P, uint32_t** _col_P, double** _val_P, uint32_t N) {

  // Get sparse matrix
  uint32_t* row_P = *_row_P;
  uint32_t* col_P = *_col_P;
  double* val_P = *_val_P;

  // Count number of elements and row counts of symmetric matrix
  int* row_counts = (int*) calloc(N, sizeof(int));
  if(row_counts == nullptr) { printf("Memory allocation failed!\n"); exit(1); }
  for(uint32_t n = 0; n < N; n++) {
    for(uint32_t i = row_P[n]; i < row_P[n + 1]; i++) {

      // Check whether element (col_P[i], n) is present
      bool present = false;
      for(uint32_t m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
	if(col_P[m] == n) {
	  present = true;
	}
      }
      if(present) row_counts[n]++;
      else {
	row_counts[n]++;
	row_counts[col_P[i]]++;
      }
    }
  }
  uint32_t no_elem = 0;
  for(uint32_t n = 0; n < N; n++) {
    no_elem += row_counts[n];
  }

  // Allocate memory for symmetrized matrix
  uint32_t* sym_row_P = (uint32_t*) malloc((N + 1) * sizeof(uint32_t));
  uint32_t* sym_col_P = (uint32_t*) malloc(no_elem * sizeof(uint32_t));
  double* sym_val_P = (double*) malloc(no_elem * sizeof(double));
  if(sym_row_P == nullptr || sym_col_P == nullptr || sym_val_P == nullptr) { printf("Memory allocation failed!\n"); exit(1); }

  // Construct new row indices for symmetric matrix
  sym_row_P[0] = 0;
  for(uint32_t n = 0; n < N; n++) sym_row_P[n + 1] = sym_row_P[n] + (uint32_t) row_counts[n];

  // Fill the result matrix
  int* offset = (int*) calloc(N, sizeof(int));
  if(offset == nullptr) { printf("Memory allocation failed!\n"); exit(1); }
  for(uint32_t n = 0; n < N; n++) {
    for(uint32_t i = row_P[n]; i < row_P[n + 1]; i++) {
      // considering element(n, col_P[i])

      // Check whether element (col_P[i], n) is present
      bool present = false;
      for(uint32_t m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
	if(col_P[m] == n) {
	  present = true;
	  if(n <= col_P[i]) {
	    // make sure we do not add elements twice
	    sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
	    sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
	    sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
	    sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
	  }
	}
      }

      // If (col_P[i], n) is not present, there is no addition involved
      if(!present) {
	sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
	sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
	sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
	sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
      }

      // Update offsets
      if(!present || (present && n <= col_P[i])) {
	offset[n]++;
	if(col_P[i] != n) offset[col_P[i]]++;
      }
    }
  }

  // Divide the result by two
  for(uint32_t i = 0; i < no_elem; i++) sym_val_P[i] /= 2.0;

  // Return symmetrized matrices
  free(*_row_P); *_row_P = sym_row_P;
  free(*_col_P); *_col_P = sym_col_P;
  free(*_val_P); *_val_P = sym_val_P;

  // Free up some memery
  free(offset); offset = nullptr;
  free(row_counts); row_counts  = nullptr;
}

// Compute squared Euclidean distance matrix
void TSNE::computeSquaredEuclideanDistance(double* X, uint32_t N, uint32_t D, double* DD) {
  const double* XnD = X;
  for(uint32_t n = 0; n < N; ++n, XnD += D) {
    const double* XmD = XnD + D;
    double* curr_elem = &DD[n*N + n];
    *curr_elem = 0.0;
    double* curr_elem_sym = curr_elem + N;
    for(uint32_t m = n + 1; m < N; ++m, XmD+=D, curr_elem_sym+=N) {
      *(++curr_elem) = 0.0;
      for(uint32_t d = 0; d < D; ++d) {
	*curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]);
      }
      *curr_elem_sym = *curr_elem;
    }
  }
}


// Makes data zero-mean
void TSNE::zeroMean(double* X, uint32_t N, uint32_t D) {

  // Compute data mean
  double* mean = (double*) calloc(D, sizeof(double));
  if(mean == nullptr) { printf("Memory allocation failed!\n"); exit(1); }
  uint32_t nD = 0;
  for(uint32_t n = 0; n < N; n++) {
    for(uint32_t d = 0; d < D; d++) {
      mean[d] += X[nD + d];
    }
    nD += D;
  }
  for(uint32_t d = 0; d < D; d++) {
    mean[d] /= (double) N;
  }

  // Subtract data mean
  nD = 0;
  for(uint32_t n = 0; n < N; n++) {
    for(uint32_t d = 0; d < D; d++) {
      X[nD + d] -= mean[d];
    }
    nD += D;
  }
  free(mean); mean = nullptr;
}


// Generates a Gaussian random number
double TSNE::randn() {
  double x, y, radius;
  do {
    x = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
    y = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
    radius = (x * x) + (y * y);
  } while((radius >= 1.0) || (radius == 0.0));
  radius = sqrt(-2 * log(radius) / radius);
  x *= radius;
  y *= radius;
  return x;
}

// Function that loads data from a t-SNE file
// Note: this function does a malloc that should be freed elsewhere
bool TSNE::load_data(double** data, int* n, int* d, int* no_dims, double* theta, double* perplexity, int* rand_seed, int* max_iter) {

  // Open file, read first 2 integers, allocate memory, and read the data
  FILE *h;
  if((h = fopen("data.dat", "r+b")) == nullptr) {
    printf("Error: could not open data file.\n");
    return false;
  }
  fread(n, sizeof(uint32_t), 1, h);	// number of datapoints
  fread(d, sizeof(uint32_t), 1, h);  // original dimensionality
  fread(theta, sizeof(double), 1, h);  // gradient accuracy
  fread(perplexity, sizeof(double), 1, h);  // perplexity
  fread(no_dims, sizeof(uint32_t), 1, h);  // output dimensionality
  fread(max_iter, sizeof(uint32_t), 1,h);  // maximum number of iterations
  *data = (double*) malloc(*d * *n * sizeof(double));
  if(*data == nullptr) {
    printf("Memory allocation failed!\n"); exit(1);
  }
  fread(*data, sizeof(double), *n * *d, h);  // the data
  if(!feof(h)) {
    fread(rand_seed, sizeof(int), 1, h); // random seed
  }
  fclose(h);
  printf("Read the %i x %i data matrix successfully!\n", *n, *d);
  return true;
}

// Function that saves map to a t-SNE file
void TSNE::save_data(double* data, int* landmarks, double* costs, uint32_t n, uint32_t d) {

  // Open file, write first 2 integers and then the data
  FILE *h;
  if((h = fopen("result.dat", "w+b")) == nullptr) {
    printf("Error: could not open data file.\n");
    return;
  }
  fwrite(&n, sizeof(int), 1, h);
  fwrite(&d, sizeof(int), 1, h);
  fwrite(data, sizeof(double), n * d, h);
  fwrite(landmarks, sizeof(int), n, h);
  fwrite(costs, sizeof(double), n, h);
  fclose(h);
  printf("Wrote the %i x %i data matrix successfully!\n", n, d);
}

/*

// Function that runs the Barnes-Hut implementation of t-SNE
int main() {

// Define some variables
uint32_t origN, N, D, no_dims, max_iter, *landmarks;
double perc_landmarks;
double perplexity, theta, *data;
int rand_seed = -1;
TSNE* tsne = new TSNE();

// Read the parameters and the dataset
if(tsne->load_data(&data, &origN, &D, &no_dims, &theta, &perplexity, &rand_seed, &max_iter)) {

// Make dummy landmarks
N = origN;
int* landmarks = (int*) malloc(N * sizeof(int));
if(landmarks == nullptr) {
printf("Memory allocation failed!\n"); exit(1);
}
for(uint32_t n = 0; n < N; n++) {
landmarks[n] = n;
}

// Now fire up the SNE implementation
double* Y = (double*) malloc(N * no_dims * sizeof(double));
double* costs = (double*) calloc(N, sizeof(double));
if(Y == nullptr || costs == nullptr) {
printf("Memory allocation failed!\n"); exit(1);
}
tsne->run(data, N, D, Y, no_dims, perplexity, theta, rand_seed, false, max_iter);

// Save the results
tsne->save_data(Y, landmarks, costs, N, no_dims);

// Clean up the memory
free(data); data = nullptr;
free(Y); Y = nullptr;
free(costs); costs = nullptr;
free(landmarks); landmarks = nullptr;
}
delete(tsne);
}
*/

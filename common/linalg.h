#ifndef _linalg_h_included_
#define _linalg_h_included_

#include <boost/config.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <jsc/util/log.hpp>

extern "C" {
	extern void dgetrf_(int *,int *,double *,int *,double *,int *);
	extern void dgetri_(int *,double *,int *,double *,double *,int *,int *);
	extern void dpotrf_(char *,int *,double *,int *,int *);
	extern void dpotri_(char *,int *,double *,int *,int *);
};

using namespace boost;
using namespace boost::numeric;
using namespace jsc::util;

class linalg {
	public:
		static void invert_matrix(const ublas::matrix<double> & A,
				ublas::matrix<double> & inv) {
			invert_matrix_gsl_lu(A, inv);
			// invert_matrix_lapack_lu(A, inv);
			// invert_matrix_lapack_cholesky(A, inv);
			unsigned long K = A.size1() + 1;
			ublas::matrix<double> rst(K-1, K-1);
			ublas::matrix<double> ident = ublas::identity_matrix<double>(K-1);
			axpy_prod(A, inv, rst);
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					if (abs(rst(p,q) - ident(p,q)) > 1E-5) {
						L_(debug) << "A*(invA)[" << p << "," << q
							<< "] = " << rst(p,q);
					}
				}
			}
		}
	private:
		static void invert_matrix_gsl_lu(const ublas::matrix<double> & A,
				ublas::matrix<double> & inv) {
			unsigned long K = A.size1() + 1;

			gsl_matrix* m = gsl_matrix_alloc(K - 1, K - 1);
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					gsl_matrix_set(m, p, q, A(p,q));
				}
			}

			gsl_matrix* inverse = gsl_matrix_alloc(K - 1, K - 1);
			gsl_permutation*perm = gsl_permutation_alloc(K-1);
			int s=0;
			gsl_linalg_LU_decomp(m, perm, &s);
			gsl_linalg_LU_invert(m, perm, inverse);
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					inv(p,q) = gsl_matrix_get(inverse, p, q);
				}
			}
			gsl_permutation_free(perm);
			gsl_matrix_free(inverse);
			gsl_matrix_free(m);
		}

		static void invert_matrix_lapack_lu(const ublas::matrix<double> & A,
				ublas::matrix<double> & inv) {
			unsigned long K = A.size1() + 1;

			double * m = (double *)malloc((K - 1)*(K - 1)*sizeof(double));
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					m[p + q*(K-1)] = A(p,q);
				}
			}
			double * ipiv = (double *)malloc((K - 1)*sizeof(double));
			int N = K - 1;
			int lda = N;
			int info;
			dgetrf_(&N, &N, &m[0], &lda, &ipiv[0], &info);
			assert(info == 0);

			double * work = (double *)malloc((K - 1)*sizeof(double));
			int lwork = N;
			dgetri_(&N, &m[0], &lda, &ipiv[0], &work[0], &lwork, &info);
			assert(info == 0);

			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					inv(p,q) = m[p + q*(K-1)];
				}
			}

			free(m);
			free(ipiv);
			free(work);
		}

		static void invert_matrix_lapack_cholesky(const ublas::matrix<double> & A,
				ublas::matrix<double> & inv) {
			unsigned long K = A.size1() + 1;

			double * m = (double *)malloc((K - 1)*(K - 1)*sizeof(double));
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					m[p + q*(K-1)] = A(p,q);
				}
			}
			int N = K - 1;
			int lda = N;
			int info;
			char UPLO = 'U';
			dpotrf_(&UPLO, &N, &m[0], &lda, &info);
			assert(info == 0);
			dpotrf_(&UPLO, &N, &m[0], &lda, &info);
			assert(info == 0);

			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					inv(p,q) = m[p + q*(K-1)];
				}
			}

			free(m);
		}
};

#endif

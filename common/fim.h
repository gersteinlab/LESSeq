#ifndef _fim_h_included_
#define _fim_h_included_

#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <vector>
#include <sstream>
#include <string>
#include <boost/config.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <boost/numeric/ublas/io.hpp>

#include <jsc/bioinfo/gene_anno.hpp>
#include <jsc/util/log.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "splicing_graph.h"
#include "read.h"
#include "linalg.h"
#include "accessible_read_starts.h"

using namespace std;
using namespace boost;
using namespace boost::numeric;
using namespace jsc::bioinfo;
using namespace jsc::util;

class FIM {
	public:
		FIM () {
			ofim_call_count = 0;
		}

		double get_ofim_call_count() const {
			return ofim_call_count;
		};

		static double estimate_mle_variance(
				ublas::matrix<double> const & I) {
			return estimate_mle_variance_by_diag(I);
		}

		static double estimate_mle_variance_by_diag(
				ublas::matrix<double> const & I) {
			unsigned long K = I.size1() + 1;
			double sum = 0;
			for (unsigned long p = 0; p < K - 1; ++p) {
				sum += 1.0 / I(p,p);
			}
			return sum;
		}

		static double estimate_mle_variance_by_inv(
				ublas::matrix<double> const & I) {
			assert(I.size1() == I.size2());
			L_(debug) << I;
			ublas::matrix<double> inv =
				ublas::zero_matrix<double>(I.size1(), I.size1());
			linalg::invert_matrix(I, inv);

			L_(debug) << inv;

			unsigned long K = I.size1() + 1;
			double sum = 0;
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					sum += inv(p,q);
				}
				sum += inv(p,p);
			}

			return sum;
		}

		static void print_fim_for_R(ostream & os, ublas::matrix<double> I) {
			unsigned long K = I.size1() + 1;
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					os << I(p,q) << "\t";
				}
				os << endl;
			}
		}

		static void print_fim_diag_for_R(ostream & os, ublas::matrix<double> I) {
			unsigned long K = I.size1() + 1;
			os << "(";
			for (unsigned long p = 0; p < K - 1; ++p) {
				os << I(p,p) << ",";
			}
			os << ")";
			os << endl;
		}

		ublas::matrix<double> bruteforce_fim(shared_ptr<Read> readp,
				vector<double> const & iso_probs) {
			L_(debug) << "FIM::bruteforce_fim";
			reset_ofim_call_count();
			shared_ptr<AccessibleReadStarts> ars = readp->get_ARS();
			unsigned long K = ars->get_num_isoforms();
			ublas::matrix<double> I = ublas::zero_matrix<double>(K-1, K-1);

			for (unsigned long k = 0; k < K; ++k) {
				if (iso_probs[k] == 0) {
					continue;
				}
				ublas::matrix<int> sign;
				L_(debug) << "Possible Isoform #" << k
					<< ": Generating reads with starts from 0 to "
					<< ars->get_iso_ARS_total_length(k)
					<< " in ARS";
				for (unsigned long a = 0;
						a < ars->get_iso_ARS_total_length(k);
						++a) {
					// generate a read starting from a in ARS
					readp->generate_read(k, ars->ARStart2IsoStart(k, a));
					double G = readp->prob_generated_by_iso(k);
					double log_scaler = log(G) + log(iso_probs[k]);
					ublas::matrix<double> m = ofim(readp, iso_probs, sign);
					for (unsigned long p = 0; p < K - 1; ++p) {
						for (unsigned long q = 0; q < K - 1; ++q) {
							if (sign(p,q) > 0) {
								m(p,q) = exp(m(p,q) + log_scaler);
							} else if (sign(p,q) < 0) {
								m(p,q) = - exp(m(p,q) + log_scaler);
							} else {
								// m(p,q) must be zero, no need to perform scaling
							}
							assert(!isinf(m(p,q)) && !isnan(m(p,q)));
						}
					}
					I += m;
				}
			}
			return I;
		}

		ublas::matrix<double> fast_fim(shared_ptr<Read> readp,
				vector<double> const & iso_probs) {
			L_(debug) << "FIM::fast_fim";
			reset_ofim_call_count();
			shared_ptr<AccessibleReadStarts> ars = readp->get_ARS();
			unsigned long K = ars->get_num_isoforms();
			ublas::matrix<double> I = ublas::zero_matrix<double>(K-1, K-1);
			for (unsigned long k = 0; k < K; ++k) {
				if (iso_probs[k] == 0) {
					continue;
				}
				ublas::matrix<int> sign;
				L_(debug) << "Possible Isoform #" << k
					<< ": Generating reads with starts from 0 to "
					<< ars->get_iso_ARS_total_length(k)
					<< " in ARS";

				unsigned long a = 0;
				while (a < ars->get_iso_ARS_total_length(k)) {
					// generate a read starting from a
					readp->generate_read(k, ars->ARStart2IsoStart(k, a));

					unsigned long N_overlapping_exons = readp->exon_indices.size();
					unsigned long N_eq_samples = min(
							ars->exon_lengths[readp->exon_indices[0]]
							- readp->start_at_first_exon,
							ars->exon_lengths[readp->exon_indices[N_overlapping_exons - 1]]
							- readp->end_at_last_exon + 1);

					while ( (a + N_eq_samples - 1 >= ars->get_iso_ARS_total_length(k)) || (ars->ARStart2IsoStart(k, a + N_eq_samples - 1) != (ars->ARStart2IsoStart(k, a) + N_eq_samples - 1)) ){
						N_eq_samples--;
					}

					double G = readp->prob_generated_by_iso(k);
					ublas::matrix<double> m = ofim(readp, iso_probs, sign);
					bool zero_scaler = ((G == 0) || (iso_probs[k] == 0));
					double log_scaler = log(N_eq_samples)
						+ log(G)
						+ log(iso_probs[k]);
					for (unsigned long p = 0; p < K - 1; ++p) {
						for (unsigned long q = 0; q < K - 1; ++q) {
							if (zero_scaler) {
								m(p,q) = 0;
							} else if (sign(p,q) > 0) {
								m(p,q) = exp(m(p,q) + log_scaler);
							} else if (sign(p,q) < 0) {
								m(p,q) = - exp(m(p,q) + log_scaler);
							} else {
								// m(p,q) must be zero, no need to perform scaling
							}
							assert(!isinf(m(p,q)) && !isnan(m(p,q)));
						}
					}
					I += m;
					a += N_eq_samples;
				}
			}
			return I;
		}

		ublas::matrix<double> faster_fim(shared_ptr<Read> readp,
				vector<double> const & iso_probs) {
			L_(debug) << "FIM::faster_fim";
			reset_ofim_call_count();
			shared_ptr<AccessibleReadStarts> ars = readp->get_ARS();
			unsigned long K = ars->get_num_isoforms();
			ublas::matrix<double> I = ublas::zero_matrix<double>(K-1, K-1);
			vector<interval_list<long > > covered_sample_starts;
			for (unsigned long k = 0; k < K; ++k) {
				interval_list<long> il;
				covered_sample_starts.push_back(il);
			}
			for (unsigned long k = 0; k < K; ++k) {
				if (iso_probs[k] == 0) {
					continue;
				}

				unsigned long a = covered_sample_starts[k].find_min_uncovered(0);
				L_(debug) << "Possible Isoform #" << k
					<< ": Generating reads with starts from " << a << " to "
					<< ars->get_iso_ARS_total_length(k)
					<< " in ARS";
				while (a < ars->get_iso_ARS_total_length(k)) {
					// generate a read starting from a
					readp->generate_read(k, ars->ARStart2IsoStart(k, a));

					unsigned long N_overlapping_exons = readp->exon_indices.size();
					unsigned long N_eq_samples = min(
							ars->exon_lengths[readp->exon_indices[0]]
							- readp->start_at_first_exon,
							ars->exon_lengths[readp->exon_indices[N_overlapping_exons - 1]]
							- readp->end_at_last_exon + 1);
					while ( (a + N_eq_samples - 1 >= ars->get_iso_ARS_total_length(k)) || (ars->ARStart2IsoStart(k, a + N_eq_samples - 1) != (ars->ARStart2IsoStart(k, a) + N_eq_samples - 1)) ){
						N_eq_samples--;
					}

					ublas::matrix<int> sign;
					ublas::matrix<double> m = ofim(readp, iso_probs, sign);

					double G = readp->prob_generated_by_iso(k);
					double scaler = G * N_eq_samples * iso_probs[k];

					if (covered_sample_starts[k].compute_overlap(a, a + N_eq_samples)
							> 0) {
						L_(debug2) << "Cov for iso # " << k
							<< ": "	<< covered_sample_starts[k];
						L_(error) << "Overlap > 0!";
					}
					covered_sample_starts[k].add_interval(a, a + N_eq_samples);
					L_(debug2) << "Adding [" << a
						<< "," << (a + N_eq_samples)
						<< ") to possible isoform #" << k;
					for (unsigned long j = 0; j < K; ++j) {
						if (j == k) continue;
						long first_matching_idx =
							is_connected_exons_compatible_with_isoform(
									readp->exon_indices,
									ars->iso_exon_indices[j]);
						if (first_matching_idx >= 0) {
							double G_j = readp->prob_generated_by_iso(j);
							scaler += G_j * N_eq_samples * iso_probs[j];
							unsigned long a2 = readp->start_at_first_exon
								+ ars->iso_exon_total_lengths[j][first_matching_idx]
								- ars->exon_lengths[readp->exon_indices[0]];
							a2 = ars->IsoStart2ARStart(j, a2);

							if (covered_sample_starts[j].compute_overlap(a2, a2 + N_eq_samples)
									> 0) {
								L_(debug2) << "Cov for iso # " << j
									<< ": "	<< covered_sample_starts[j];
								L_(error) << "Overlap > 0!";
							}
							covered_sample_starts[j].add_interval(a2, a2 + N_eq_samples);
							L_(debug2) << "Adding [" << a2
								<< "," << (a2 + N_eq_samples)
								<< ") to possible isoform #" << j;
						}
					}

					double log_scaler = log(scaler);
					for (unsigned long p = 0; p < K - 1; ++p) {
						for (unsigned long q = 0; q < K - 1; ++q) {
							if (sign(p,q) > 0) {
								m(p,q) = exp(m(p,q) + log_scaler);
							} else if (sign(p,q) < 0) {
								m(p,q) = - exp(m(p,q) + log_scaler);
							} else {
								// m(p,q) must be zero, no need to perform scaling
							}
							assert(!isinf(m(p,q)) && !isnan(m(p,q)));
						}
					}
					I += m;
					a = covered_sample_starts[k].find_min_uncovered(0);
				}
			}
			return I;
		}
	protected:
		void reset_ofim_call_count() {
			ofim_call_count = 0;
		};

		// \mathfrak{I}_{s}^{(m)}(\Theta)
		ublas::matrix<double> ofim(shared_ptr<const Read> readp,
				vector<double> const & iso_probs,
				ublas::matrix<int> & sign) {
			++ofim_call_count;
			shared_ptr<AccessibleReadStarts> ars = readp->get_ARS();
			unsigned long K = ars->get_num_isoforms();
			ublas::matrix<double> m = ublas::zero_matrix<double>(K-1,K-1);
			sign = ublas::zero_matrix<int>(K-1, K-1);
			vector<double> v_delta_G(K, 0);
			for (unsigned long k = 0; k < K; ++k) {
				if (readp->is_compatible_with_iso(k)) {
					v_delta_G[k] = readp->prob_generated_by_iso(k);
				}
			}

			double sum = 0;
			for (unsigned long k = 0; k < K; ++k) {
				if (v_delta_G[k] > 0 && iso_probs[k] > 0) {
					sum += iso_probs[k] * v_delta_G[k];
				}
			}

			assert(sum > 0);

			double log_prod = log(sum) + log(sum);
			for (unsigned long p = 0; p < K - 1; ++p) {
				for (unsigned long q = 0; q < K - 1; ++q) {
					m(p,q) = (v_delta_G[p] - v_delta_G[K - 1])*
						(v_delta_G[q] - v_delta_G[K - 1]);
					if (m(p,q) != 0) {
						sign(p,q) = m(p,q) > 0 ? 1 : -1;
						if (m(p,q) < 0) {
							m(p,q) = - m(p,q);
						}
						m(p,q) = log(m(p,q)) - log_prod;
					}
				}
			}

			if (readp->exon_indices[0] == 0
					&& readp->start_at_first_exon == 2730) {
				L_(info) << *readp << " " << m(0,0) << " " << sign(0,0);
			}

			return m;
		}

	private:
		double ofim_call_count;
};

#endif

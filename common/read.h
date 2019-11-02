#ifndef _read_h_included_
#define _read_h_included_

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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <jsc/bioinfo/gene_anno.hpp>
#include <jsc/util/log.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "splicing_graph.h"
#include "accessible_read_starts.h"

using namespace std;
using namespace boost;
using namespace boost::numeric;
using namespace jsc::bioinfo;
using namespace jsc::util;

using boost::shared_ptr;

// returns the index (in the isoform) of the first matching exon in the isoform
long is_connected_exons_compatible_with_isoform(
		vector<unsigned long> const & exon_indices,
		vector<unsigned long> const & iso_exon_indices) {
	vector<unsigned long>::const_iterator itr, itr2;
	itr = exon_indices.begin();
	itr2 = iso_exon_indices.begin();
	bool is_compatible = false;
	bool found_first_match = false;
	long first_match_idx = 0;
	while (itr != exon_indices.end()) {
		if (itr2 == iso_exon_indices.end()) {
			is_compatible = false;
			break;
		}
		if (*itr != *itr2) {
			if (!found_first_match) {
				itr2++;
				first_match_idx++;
			} else {
				is_compatible = false;
				break;
			}
		} else {
			if (!found_first_match) {
				found_first_match = true;
				is_compatible = true;
				itr++;
				itr2++;
			} else {
				itr++;
				itr2++;
			}
		}
	}
	return (is_compatible ? first_match_idx : -1);
}

class ReadType {
	public:
		static const int LONG_READ = 0;
		static const int MEDIUM_READ = 1;
		static const int SHORT_READ = 2;
		static const int SHORT_CONTIG = 3;
		static const int CONTIG = 4;
		static const int LONG_READ_PE = 5;
		static const int MEDIUM_READ_PE = 6;
		static const int SHORT_READ_PE = 7;
		static const string TYPES[];
};
const string ReadType::TYPES[] = {"LONG_READ", "MEDIUM_READ", "SHORT_READ", "SHORT_CONTIG", "CONTIG", "LONG_READ_PE", "MEDIUM_READ_PE", "SHORT_READ_PE"};

class Read {
	public:
		string read_type;
		vector<unsigned long> exon_indices;
		unsigned long start_at_first_exon;
		unsigned long end_at_last_exon;

		Read(shared_ptr<AccessibleReadStarts> const & ars
				= shared_ptr<AccessibleReadStarts>())
			: ars(ars) {
			read_length = 0;
			expected_read_length = 0;
			read_type = "General read";
		}

		unsigned long get_num_isoforms() const {
			return ars->get_num_isoforms();
		}

		unsigned long get_read_length() const {
			return read_length;
		}

		shared_ptr<AccessibleReadStarts> get_ARS() const {
			return ars;
		}

		virtual unsigned long get_read_span() const {
			return read_length;
		}

		unsigned long get_expected_read_length() const {
			return expected_read_length;
		}

		virtual unsigned long get_expected_read_span() const {
			return expected_read_length;
		}

		virtual bool is_compatible_with_iso(
				unsigned long const & iso_idx) const = 0;

		virtual void build(
				interval_list<long> const & read_il,
				SetExonp const & exonps,
				interval_list<long> const & read_il2 = interval_list<long>()) = 0;

		/// Generates a read from an isoform.
		//
		/// The start position of the read in the isoform is either
		/// determined by \a fixed_read_start or chosen randomly
		/// according to the actual type of the read.
		//
		/// \param[in] iso_idx							index of the isoform from which
		///																	the read will be generated.
		///																	exons in the corresponding gene.
		/// \param[in] fixed_read_start			When >= 0, specifies the start
		///																	position (in the isoform) of the
		///																	read to be generated.
		///	\returns	The actual start position (in the isoform)
		///						of the generated read.
		virtual unsigned long generate_read(
				unsigned long const & iso_idx,
				long fixed_read_start = -1) = 0;

		virtual double prob_generated_by_iso(
				unsigned long const & iso_idx) const = 0;

		virtual string to_string() const = 0;

		virtual ~Read() {}

	protected:
		shared_ptr<AccessibleReadStarts> ars;
		unsigned long read_length;
		unsigned long expected_read_length;
};

ostream& operator << (ostream& os, Read const & r) {
	os << r.to_string();
	return os;
};

typedef shared_ptr<Read> Read_ptr;

class Read_stats {
	public:
		static const unsigned long medium_single_expected_read_length = 250;
		static const double medium_cost_per_bp = 7E-5;
		static const unsigned long short_single_expected_read_length = 30;
		static const unsigned long short_paired_expected_insert_size = 300;
		static const double short_cost_per_bp = 7E-6;
};

class Read_single : public Read {
	public:
		Read_single(shared_ptr<AccessibleReadStarts> ars
				= shared_ptr<AccessibleReadStarts>())
			: Read(ars)	{
			rng = NULL;
			read_type = "General single read";
		}

		virtual bool is_compatible_with_iso(
				unsigned long const & iso_idx) const {
			return (is_connected_exons_compatible_with_isoform(
						exon_indices, ars->iso_exon_indices[iso_idx]) >= 0);
		}

		virtual void build(
				interval_list<long> const & read_il,
				SetExonp const & exonps,
				interval_list<long> const & read_il2 = interval_list<long>()) {
			L_(debug2) << "Read::build";
			unsigned long matching_length = 0;
			set<unsigned long> index_set;
			exon_indices.clear();
			start_at_first_exon = end_at_last_exon = 0;
			bool found_start = false;
			unsigned long il_idx = 0;
			SetExonp::const_iterator itr = exonps.begin();
			long cur_start, cur_end, cur_exon_start;
			cur_exon_start = 0;
			L_(debug2) << "Read: " << read_il;
			L_(debug2) << "Exons: " << exonps;
			while (il_idx < read_il.get_num_intervals()) {
				cur_start = read_il.get_starts()[il_idx];
				cur_end = read_il.get_ends()[il_idx];
				L_(debug2) << "Current interval: [" << cur_start << ","
					<< cur_end << ")";
				// find all the exons that overlap with/contain the interval
				while (itr != exonps.end() && (*itr)->start < cur_end) {
					L_(debug2) << "Current exon: [" << (*itr)->start << ","
						<< (*itr)->end << ")";
					cur_exon_start = max(cur_exon_start, (*itr)->start);
					if (cur_start >= cur_exon_start && cur_start < (*itr)->end) {
						if (!found_start) {
							L_(debug2) << "Found start";
							found_start = true;
							start_at_first_exon = cur_start - (*itr)->start;
						} else if (cur_start > cur_exon_start) {
							break;
						} 

						L_(debug2) << "Add exon #" << distance(exonps.begin(), itr);
						index_set.insert(distance(exonps.begin(), itr));
						cur_exon_start = min((*itr)->end, cur_end);
						matching_length += cur_exon_start - cur_start;

						if (cur_end < (*itr)->end) {
							cur_start = cur_end;
							end_at_last_exon = cur_end - (*itr)->start;
							++il_idx;
							break;
						} else if (cur_end == (*itr)->end) {
							cur_start = cur_end;
							end_at_last_exon = cur_end - (*itr)->start;
							++il_idx;
						} else {
							cur_start = (*itr)->end;
							end_at_last_exon = (*itr)->end - (*itr)->start;
						}
					} else if (cur_exon_start > (*itr)->start
							&& cur_exon_start < (*itr)->end) {
						break;
					}
					itr++;
				}
				if (cur_start == cur_end) {
					continue;
				} else {
					break;
				}
			}

			BOOST_FOREACH (unsigned long i, index_set) {
				exon_indices.push_back(i);
			}
			read_length = matching_length;
		}

		virtual unsigned long generate_read(
				unsigned long const & iso_idx,
				long fixed_read_start = -1) {
			assert(rng != NULL || fixed_read_start >= 0);
			vector<unsigned long> const & iso_exon_indices =
				ars->iso_exon_indices[iso_idx];
			vector<unsigned long> const & iso_exon_total_lengths =
				ars->iso_exon_total_lengths[iso_idx];

			exon_indices.clear();
			unsigned long read_start;
			if (fixed_read_start < 0) {	// randomly generate a read start
				read_start = (unsigned long)(
						gsl_ran_flat(rng, 0,
							ars->get_iso_ARS_total_length(iso_idx)));
				read_start = ars->ARStart2IsoStart(iso_idx, read_start);
			} else {
				read_start = fixed_read_start;
			}

			unsigned long read_end = read_start + get_expected_read_length();
			read_length = read_end - read_start;
			L_(debug2) << "Generating read [" << read_start
				<< ", " << read_end << ") length=" << read_length << "bp";
			vector<unsigned long>::const_iterator itr =
				upper_bound(iso_exon_total_lengths.begin(),
						iso_exon_total_lengths.end(),
						read_start);
			assert(itr != iso_exon_total_lengths.end());
			unsigned long start_exon_idx = distance(
					iso_exon_total_lengths.begin(), itr);
			itr = lower_bound(iso_exon_total_lengths.begin(),
					iso_exon_total_lengths.end(),
					read_end);
			assert(itr != iso_exon_total_lengths.end());
			unsigned long end_exon_idx = distance(
					iso_exon_total_lengths.begin(), itr);
			L_(debug2) << "Start exon: "
				<< iso_exon_indices[start_exon_idx]
				<< ", end exon: "
				<< iso_exon_indices[end_exon_idx];
			for (unsigned long i = start_exon_idx;
					i <= end_exon_idx; ++i) {
				exon_indices.push_back(iso_exon_indices[i]);
			}
			start_at_first_exon = read_start -
				(start_exon_idx == 0 ?
				 0 : iso_exon_total_lengths[start_exon_idx - 1]);
			end_at_last_exon = read_end -
				(end_exon_idx == 0 ?
				 0 : iso_exon_total_lengths[end_exon_idx - 1]);

			return read_start;
		}

		virtual double prob_generated_by_iso(
				unsigned long const & iso_idx) const {
			double num_diff_reads = ars->get_iso_ARS_total_length(iso_idx);
			if (num_diff_reads <= 0) {
				L_(warning) << read_type << ": supported read length > isoform length" << endl;
				return 0;
			}
			double prob = (double)1.0 / num_diff_reads;
			return prob;
		}

		virtual string to_string() const {
			std::ostringstream os;
			os << "[" << read_type << "] "
				<< "(Exons) ";
			BOOST_FOREACH(unsigned long const & idx, exon_indices) {
				os << idx << ",";
			}
			os << " (Start at first exon) " << start_at_first_exon << " ";
			os << "(End at last exon) " << end_at_last_exon;
			return os.str();
		}
	protected:
		gsl_rng * rng;
};

class Read_medium_single : public Read_single {
	public:
		Read_medium_single(shared_ptr<AccessibleReadStarts> const & ars
				= shared_ptr<AccessibleReadStarts>(),
				gsl_rng * r = NULL,
				unsigned long expected_length =
				Read_stats::medium_single_expected_read_length)
			: Read_single(ars) {
			rng = r;
			expected_read_length = expected_length;
			read_type = "Medium single read";
		}
};

class Read_short_single : public Read_single {
	public:
		Read_short_single(shared_ptr<AccessibleReadStarts> const & ars
				= shared_ptr<AccessibleReadStarts>(),
				gsl_rng * r = NULL,
				unsigned long expected_length =
				Read_stats::short_single_expected_read_length)
			: Read_single(ars) {
			rng = r;
			expected_read_length = expected_length;
			read_type = "Short single read";
		}
};

class Read_paired: public Read {
	public:
		Read_single end1;
		Read_single end2;
		Read_single span;

		Read_paired(shared_ptr<AccessibleReadStarts> ars
				= shared_ptr<AccessibleReadStarts>(),
				double tolerance = 0)
			: Read(ars), tolerance(tolerance) {
			insert_size = expected_insert_size = 0;
			read_type = "General paired-end read";
			end1 = Read_single(ars);
			end2 = Read_single(ars);
			span = Read_single(ars);
		}

		unsigned long get_expected_insert_size() const {
			return expected_insert_size;
		}

		virtual unsigned long get_read_span() const {
			return read_length + insert_size;
		}

		virtual unsigned long get_expected_read_span() const {
			return (end1.get_expected_read_length()
				+ end2.get_expected_read_length()
				+ expected_insert_size);
		}

		virtual bool is_compatible_with_iso(
				unsigned long const & iso_idx) const {
			assert(end1.exon_indices.size() > 0
					&& end2.exon_indices.size() > 0);
			if (end1.is_compatible_with_iso(iso_idx)
					&& end2.is_compatible_with_iso(iso_idx)) {
				if (expected_insert_size == 0) {
					return true;
				}
				vector<unsigned long> const & iso_exon_indices =
					ars->iso_exon_indices[iso_idx];
				vector<unsigned long> const & exon_total_lengths =
					ars->iso_exon_total_lengths[iso_idx];

				vector<unsigned long>::const_iterator itr;
				unsigned long end1_last_exon_idx = end1.exon_indices
					[end1.exon_indices.size() - 1];
				itr = lower_bound(iso_exon_indices.begin(),
						iso_exon_indices.end(),
						end1_last_exon_idx);
				unsigned long exon_idx = distance(
						iso_exon_indices.begin(), itr);
				unsigned long end1_end =
					(exon_idx > 0 ? exon_total_lengths[exon_idx - 1] : 0)
					+ end1.end_at_last_exon;
				unsigned long end2_first_exon_idx = end2.exon_indices[0];
				itr = lower_bound(iso_exon_indices.begin(),
						iso_exon_indices.end(),
						end2_first_exon_idx);
				exon_idx = distance(iso_exon_indices.begin(), itr);
				unsigned long end2_start =
					(exon_idx > 0 ? exon_total_lengths[exon_idx - 1] : 0)
					+ end2.start_at_first_exon;
				double gap = end2_start - end1_end;
				double diff = gap / (double)(expected_insert_size);
				if (diff <= (1.0 + tolerance + epsilon)
						&& diff >= (1.0 - tolerance - epsilon)) {
					return true;
				}
			}
			return false;
		}

		virtual double prob_generated_by_iso(
				unsigned long const & iso_idx) const {
			double num_diff_reads = ars->get_iso_ARS_total_length(iso_idx);
			if (num_diff_reads <= 0) {
				L_(warning) << read_type << ": supported paired-end read length (including insert) > isoform length" << endl;
				return 0;
			}
			double prob = (double)1.0 / num_diff_reads;
			return prob;
		}

		virtual string to_string() const {
			std::ostringstream os;
			os << "[" << read_type << "]"
				<< " End #1: " << end1.to_string()
				<< "; End #2: " << end2.to_string();
			return os.str();
		}

		virtual void build(
				interval_list<long> const & read_il,
				SetExonp const & exonps,
				interval_list<long> const & read_il2 = interval_list<long>()) {
			L_(debug2) << "Read::build";
			end1.build(read_il, exonps);
			end2.build(read_il2, exonps);
			read_length = end1.get_read_length() + end2.get_read_length();
			/// \todo compute insert size and construct span, exon_indices, etc.
		}

		virtual unsigned long generate_read(
				unsigned long const & iso_idx,
				long fixed_read_start = -1) {
			assert(rng != NULL || fixed_read_start >= 0);

			unsigned long read_start;
			if (fixed_read_start < 0) {	// randomly generate a read start
				read_start = (unsigned long)(
						gsl_ran_flat(rng, 0,
							ars->get_iso_ARS_total_length(iso_idx)));
				read_start = ars->ARStart2IsoStart(iso_idx, read_start);
			} else {
				read_start = fixed_read_start;
			}

			unsigned long end1_start =
				end1.generate_read(iso_idx, read_start);
			unsigned long end2_start =
				end1_start + end1.get_read_length()
				+ expected_insert_size;
			end2.generate_read(iso_idx, end2_start);
			read_length = end1.get_read_length() + end2.get_read_length();
			
			span.generate_read(iso_idx, read_start);
			exon_indices = span.exon_indices;
			start_at_first_exon = span.start_at_first_exon;
			end_at_last_exon = span.end_at_last_exon;

			return end1_start;
		}

	protected:
		unsigned long insert_size;
		unsigned long expected_insert_size;
		double tolerance;
		static const double epsilon = 1E-5;
		gsl_rng * rng;
};

void generate_reads(unsigned long const & num_reads,
		Isoforms const & iso,
		shared_ptr<Read> readp,
		gsl_rng * rng,
		ublas::matrix<double> & m_delta_G,
		bool known_isoforms_only = false) {
	gsl_ran_discrete_t * discrete_g = gsl_ran_discrete_preproc(
			iso.num_known_isoforms, iso.known_iso_probs);

	for (unsigned long i = 0; i < num_reads; ++i) {
		unsigned long iso_idx = gsl_ran_discrete(rng, discrete_g);
		L_(debug2) << "Selected known isoform " << iso_idx;
		if (!known_isoforms_only) {
			iso_idx = iso.known_iso_indices[iso_idx];
		}
		readp->generate_read(iso_idx);
		L_(debug2) << *readp;
		if (known_isoforms_only) {
			for (unsigned long j = 0;
					j < iso.num_known_isoforms; ++j) {
				if (readp->is_compatible_with_iso(j)) {
					m_delta_G(i,j) = readp->prob_generated_by_iso(j);
				}
			}
		} else {
			for (unsigned long j = 0;
					j < iso.num_possible_isoforms; ++j) {
				if (readp->is_compatible_with_iso(j)) {
					m_delta_G(i,j) = readp->prob_generated_by_iso(j);
				}
			}
		}
	}

	gsl_ran_discrete_free(discrete_g);
}

class Read_short_paired: public Read_paired {
	public:
		Read_short_paired(shared_ptr<AccessibleReadStarts> ars
				= shared_ptr<AccessibleReadStarts>(),
				gsl_rng * r = NULL,
				unsigned long expected_length =
				2 * Read_stats::short_single_expected_read_length,
				unsigned long expected_insert =
				Read_stats::short_paired_expected_insert_size,
				double tolerance = 0)
			: Read_paired(ars, tolerance) {
			rng = r;
			expected_read_length = expected_length;
			insert_size = expected_insert_size = expected_insert;
			read_type = "Short paired-end read";
			shared_ptr<AccessibleReadStarts> end_ars(
					new AccessibleReadStarts(
						ars->iso_exon_indices,
						ars->iso_exon_total_lengths,
						ars->exon_lengths,
						expected_length / 2));
			end1 = Read_short_single(end_ars, r, expected_length / 2);
			end2 = Read_short_single(end_ars, r, expected_length / 2);
			span = Read_medium_single(ars, r, expected_length + expected_insert);
		}
};

void compute_iso_probs_EM_step(
		unsigned long num_possible_isoforms,
		vector<ublas::matrix<double> > const & m_delta_Gs,
		vector<double> const & old_theta,
		vector<double> & new_theta) {
	for (unsigned k = 0; k < num_possible_isoforms; ++k) {
		double sum_zeta = 0;
		double num_total_reads = 0;
		BOOST_FOREACH (ublas::matrix<double> const & m_delta_G, m_delta_Gs) {
			unsigned long num_reads = m_delta_G.size1();
			num_total_reads += (double)num_reads;
			for (unsigned long i = 0; i < num_reads; ++i) {
				double sum_local_probs = 0;
				for (unsigned k2 = 0; k2 < num_possible_isoforms; ++k2) {
					sum_local_probs += old_theta[k2] * m_delta_G(i, k2);
				}
				if (sum_local_probs > 0) {
					double local_prob = old_theta[k] * m_delta_G(i, k);
					if (local_prob > 0) {
						sum_zeta += local_prob / sum_local_probs;
					}
				}
			}
		}
		new_theta[k] = sum_zeta / num_total_reads;
	}
}

double reads_log_likelihood(
		unsigned long num_possible_isoforms,
		vector<ublas::matrix<double> > const & m_delta_Gs,
		vector<double> const & theta) {
	double ll = 0;
	BOOST_FOREACH (ublas::matrix<double> const & m_delta_G, m_delta_Gs) {
		unsigned long num_reads = m_delta_G.size1();
		for (unsigned long i = 0; i < num_reads; ++i) {
			double sum_local_probs = 0;
			for (unsigned k = 0; k < num_possible_isoforms; ++k) {
				sum_local_probs += theta[k] * m_delta_G(i, k);
			}
			ll += log(sum_local_probs);
		}
	}
	return ll;
}

void compute_iso_probs_EM(
		unsigned long num_possible_isoforms,
		vector<ublas::matrix<double> > const & m_delta_Gs,
		vector<double> & theta) {
	theta.resize(num_possible_isoforms,
			(double)1.0 / (double)num_possible_isoforms);
	double ll, old_ll;
	vector<double> old_theta;
	unsigned long count = 0;
	do {
		old_theta = theta;
		old_ll = reads_log_likelihood(
				num_possible_isoforms, m_delta_Gs, old_theta);
		compute_iso_probs_EM_step(num_possible_isoforms,
				m_delta_Gs, old_theta, theta);
		ll = reads_log_likelihood(
				num_possible_isoforms, m_delta_Gs, theta);
		if (++count == 100) {
			L_(debug) << "log likelihood: " << old_ll << " -> " << ll;
			count = 0;
		}
	} while (abs(1.0 - old_ll / ll) > 1E-6);
}

void err_known_iso_probs(Isoforms const & iso,
		vector<double> const & theta, double & abs_err,
		double & sse) {
	abs_err = sse = 0;
	for (unsigned long i = 0; i < iso.num_known_isoforms; ++i) {
		abs_err += abs(iso.known_iso_probs[i]
				- theta[iso.known_iso_indices[i]]);
		sse += (iso.known_iso_probs[i] - theta[iso.known_iso_indices[i]]) *
			(iso.known_iso_probs[i] - theta[iso.known_iso_indices[i]]);
	}
}

void err_possible_iso_probs(Isoforms const & iso,
		vector<double> const & theta, double & abs_err,
		double & sse) {
	abs_err = sse = 0;
	for (unsigned long i = 0; i < iso.num_possible_isoforms; ++i) {
		abs_err += theta[i];
		sse += theta[i] * theta[i];
	}
	for (unsigned long i = 0; i < iso.num_known_isoforms; ++i) {
		abs_err -= theta[iso.known_iso_indices[i]];
		abs_err += abs(iso.known_iso_probs[i]
				- theta[iso.known_iso_indices[i]]);
		sse -= theta[iso.known_iso_indices[i]] * theta[iso.known_iso_indices[i]];
		sse += (iso.known_iso_probs[i] - theta[iso.known_iso_indices[i]]) *
			(iso.known_iso_probs[i] - theta[iso.known_iso_indices[i]]);
	}
}

#endif

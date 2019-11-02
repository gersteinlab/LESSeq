#ifndef _accessible_read_starts_h_included_
#define _accessible_read_starts_h_included_

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

#include <jsc/util/log.hpp>
#include <jsc/util/interval_list.hpp>

#include "splicing_graph.h"

class AccessibleReadStarts {
	public:
		AccessibleReadStarts() {}

		AccessibleReadStarts(
				vector<vector<unsigned long> > const & iso_exon_indices,
				vector<vector<unsigned long> > const & iso_exon_total_lengths,
				vector<unsigned long> const & exon_lengths,
				unsigned long read_length)
		: iso_exon_indices(iso_exon_indices),
		iso_exon_total_lengths(iso_exon_total_lengths),
		exon_lengths(exon_lengths),
		read_length(read_length) {
		}

		virtual ~AccessibleReadStarts() {}

		virtual void construct_ARS() {
			BOOST_FOREACH(vector<unsigned long> const & exon_indices,
					iso_exon_indices) {
				vector<interval_list<unsigned long> > ars_in_exons;
				unsigned long total_length = 0;
				vector<unsigned long> total_lengths;
				unsigned long iso_length = 0;
				unsigned long iso_total_length = 0;
				for (unsigned long i = 0;
						i < exon_indices.size(); ++i) {
					iso_total_length += exon_lengths[exon_indices[i]];
				}
				for (unsigned long i = 0;
						i < exon_indices.size(); ++i) {
					unsigned long exon_idx = exon_indices[i];
					interval_list<unsigned long> ars;
					iso_length += exon_lengths[exon_idx];
					if (iso_length + read_length > iso_total_length) {
						ars.add_interval(0, 
								(unsigned long)max((long)exon_lengths[exon_idx] + 1 -
									(long)(iso_length + read_length - iso_total_length),
									(long)0));
						total_length += (unsigned long)ars.compute_total_length();
						total_lengths.push_back(total_length);
						ars_in_exons.push_back(ars);
						L_(debug2) << ars;
						for (unsigned long j = i + 1; j < exon_indices.size(); ++j) {
							interval_list<unsigned long> empty_ars;
							total_lengths.push_back(total_length);
							ars_in_exons.push_back(empty_ars);
						}
						break;
					}
					ars.add_interval(0, exon_lengths[exon_idx]);
					total_length += (unsigned long)ars.compute_total_length();
					total_lengths.push_back(total_length);
					ars_in_exons.push_back(ars);
				}
				iso_ARS_in_exons.push_back(ars_in_exons);
				iso_ARS_total_lengths.push_back(total_lengths);
			}
		}

		unsigned long IsoStart2ARStart(
				unsigned long iso_idx,
				unsigned long iso_start) const {
			vector<unsigned long> const & ARS_total_lengths =
				iso_ARS_total_lengths[iso_idx];
			vector<unsigned long> const & exon_total_lengths =
				iso_exon_total_lengths[iso_idx];
			vector<unsigned long>::const_iterator itr =
				upper_bound(exon_total_lengths.begin(),
						exon_total_lengths.end(),
						iso_start);
			assert(itr != exon_total_lengths.end());
			unsigned long start_exon_idx = distance(
					exon_total_lengths.begin(), itr);
			L_(debug2) << "start_exon_idx: " << start_exon_idx;

			unsigned long ARS_total_length_before_start =
				start_exon_idx > 0 ? ARS_total_lengths[start_exon_idx - 1] : 0;
			unsigned long exon_total_length_before_start =
				start_exon_idx > 0 ? exon_total_lengths[start_exon_idx - 1] : 0;

			unsigned long start_at_start_exon = iso_start
				- exon_total_length_before_start;

			interval_list<unsigned long> const & start_ARS =
				iso_ARS_in_exons[iso_idx][start_exon_idx];
			unsigned long ars_total = 0;
			for (unsigned long i = 0; i < start_ARS.get_num_intervals(); ++i) {
				if (start_ARS.get_starts()[i] < start_at_start_exon) {
					ars_total += min(start_ARS.get_ends()[i], start_at_start_exon)
						- start_ARS.get_starts()[i];
				} else {
					break;
				}
			}
			unsigned long ars_start =
				ARS_total_length_before_start
				+ ars_total;
			L_(debug2) << "ars_start: " << ars_start;
			return ars_start;
		}

		unsigned long ARStart2IsoStart(
				unsigned long iso_idx,
				unsigned long start) const {
			assert(start < get_iso_ARS_total_length(iso_idx));
			vector<unsigned long> const & ARS_total_lengths =
				iso_ARS_total_lengths[iso_idx];
			vector<unsigned long> const & exon_total_lengths =
				iso_exon_total_lengths[iso_idx];
			vector<unsigned long>::const_iterator itr =
				upper_bound(ARS_total_lengths.begin(),
						ARS_total_lengths.end(),
						start);
			assert(itr != ARS_total_lengths.end());
			unsigned long start_ARS_idx = distance(
					ARS_total_lengths.begin(), itr);
			L_(debug2) << "start_ARS_idx: " << start_ARS_idx;
			unsigned long ARS_total_length_before_start =
				start_ARS_idx > 0 ? ARS_total_lengths[start_ARS_idx - 1] : 0;
			L_(debug2) << "ARS_total_length_before_start: " << ARS_total_length_before_start;
			unsigned long exon_total_length_before_start =
				start_ARS_idx > 0 ? exon_total_lengths[start_ARS_idx - 1] : 0;
			L_(debug2) << "exon_total_length_before_start: " << exon_total_length_before_start;
			unsigned long start_at_start_ARS = start - ARS_total_length_before_start;
			L_(debug2) << "start_at_start_ARS: " << start_at_start_ARS;
			interval_list<unsigned long> const & start_ARS =
				iso_ARS_in_exons[iso_idx][start_ARS_idx];
			unsigned long iso_start_at_start_exon = 0;
			unsigned long ars_total = 0;
			for (unsigned long i = 0; i < start_ARS.get_num_intervals(); ++i) {
				iso_start_at_start_exon += start_ARS.get_starts()[i]
					- (i > 0 ? start_ARS.get_ends()[i-1] : 0);
				unsigned long length = start_ARS.get_ends()[i]
					- start_ARS.get_starts()[i];
				L_(debug2) << "ars #: " << i;
				L_(debug2) << "ars_total: " << ars_total;
				if (ars_total + length <= start_at_start_ARS) {
					ars_total += length;
					iso_start_at_start_exon += length;
				} else {
					iso_start_at_start_exon += start_at_start_ARS - ars_total;
					break;
				}
			}
			L_(debug2) << "iso_start_at_start_exon: " << iso_start_at_start_exon;
			unsigned long iso_start =
				exon_total_length_before_start
				+ iso_start_at_start_exon;
			L_(debug2) << "iso_start: " << iso_start;
			return iso_start;
		}

		unsigned long get_iso_ARS_total_length(unsigned long iso_idx) const {
			assert(iso_idx < get_num_isoforms());
			unsigned long num_exons = iso_exon_indices[iso_idx].size();
			assert(num_exons > 0);
			return iso_ARS_total_lengths[iso_idx][num_exons - 1];
		}

		unsigned long get_num_isoforms() const {
			return iso_exon_indices.size();
		}

	public:
		const vector<vector<unsigned long> > iso_exon_indices;
		const vector<vector<unsigned long> > iso_exon_total_lengths;
		const vector<unsigned long> exon_lengths;
	protected:
		//! The ARS in each exon of the isoforms
		vector<vector<interval_list<unsigned long> > > iso_ARS_in_exons;
		vector<vector<unsigned long> > iso_ARS_total_lengths;
		unsigned long read_length;
};

class AccessibleShortReadStarts: public AccessibleReadStarts {
	public:
		AccessibleShortReadStarts() {}

		AccessibleShortReadStarts(
				vector<vector<unsigned long> > const & iso_exon_indices,
				vector<vector<unsigned long> > const & iso_exon_total_lengths,
				vector<unsigned long> const & exon_lengths,
				unsigned long read_length,
				unsigned long min_partial_exon_size)
		: AccessibleReadStarts(iso_exon_indices,
				iso_exon_total_lengths, exon_lengths, read_length),
			min_partial_exon_size(min_partial_exon_size) {
		}

		virtual void construct_ARS() {
			BOOST_FOREACH(vector<unsigned long> const & exon_indices,
					iso_exon_indices) {
				L_(debug2) << "constructing ars for isoform:" << exon_indices;
				vector<interval_list<unsigned long> > ars_in_exons;
				unsigned long total_length = 0;
				vector<unsigned long> total_lengths;
				unsigned long iso_length = 0;
				unsigned long iso_total_length = 0;
				for (unsigned long i = 0;
						i < exon_indices.size(); ++i) {
					iso_total_length += exon_lengths[exon_indices[i]];
				}
				for (unsigned long i = 0;
						i < exon_indices.size(); ++i) {
					unsigned long exon_idx = exon_indices[i];
					interval_list<unsigned long> ars;
					unsigned long last = max(
							(long)exon_lengths[exon_idx] - (long)read_length + 1, (long)0);

					iso_length += exon_lengths[exon_idx];
					if (iso_length + read_length > iso_total_length) {
						unsigned long last2 = max((long)exon_lengths[exon_idx] + 1 -
								(long)(iso_length + read_length - iso_total_length),
								(long)0);
						ars.add_interval(0, last2);
						total_length += (unsigned long)ars.compute_total_length();
						total_lengths.push_back(total_length);
						ars_in_exons.push_back(ars);
						for (unsigned long j = i + 1; j < exon_indices.size(); ++j) {
							interval_list<unsigned long> empty_ars;
							total_lengths.push_back(total_length);
							ars_in_exons.push_back(empty_ars);
						}
						break;
					}

					ars.add_interval(0, last);
					unsigned long partial_overlap_start = max(
							(long)exon_lengths[exon_idx] - (long)read_length
							+ (long)min_partial_exon_size, (long)0);
					if (exon_lengths[exon_idx] > min_partial_exon_size) {
						ars.add_interval(partial_overlap_start,
								exon_lengths[exon_idx] - min_partial_exon_size + 1);
					}
					total_length += (unsigned long)ars.compute_total_length();
					total_lengths.push_back(total_length);
					ars_in_exons.push_back(ars);
					L_(debug2) << "ars in exon #" << exon_idx << ": " << ars;
				}
				iso_ARS_in_exons.push_back(ars_in_exons);
				iso_ARS_total_lengths.push_back(total_lengths);
			}
		}
	protected:
		unsigned long min_partial_exon_size;
};

#endif

#ifndef _splicing_graph_h_included_
#define _splicing_graph_h_included_

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

using namespace std;
using namespace boost;
using namespace boost::numeric;
using namespace jsc::bioinfo;
using namespace jsc::util;

class Exon {
public:
	long start;
	long end;
	Exon (long a, long b)
		: start(a), end(b) {}
	bool is_inside(long a, long b) const {
		return a <= start && end <= b;
	}
};
typedef shared_ptr<Exon> ExonPtr;
typedef vector<ExonPtr> VecExonp;

ostream& operator << (ostream& os, Exon const & e)
{
	os << "[" << e.start << "," << e.end << ")";
	return os;
};

template <class T>
ostream& operator << (ostream& os, vector<T> const & v)
{
	os << "(";
	BOOST_FOREACH(T const & e, v) {
    os << e << ",";
	}
	os << ")";
	return os;
};

class ExonComp {
public:
	bool operator() (ExonPtr const & a, ExonPtr const & b) const {
		return a->start < b->start;
	}
};

typedef set<ExonPtr, ExonComp> SetExonp;

ostream& operator << (ostream& os, SetExonp const & exonps)
{
	BOOST_FOREACH (ExonPtr const & exonp, exonps){
    os << *exonp << "-";
	}
	return os;
};

class ExonSet {
public:
	// in this set, exon1->start < exon1->end < exon2->start < exon2->end...
	SetExonp exonps;

	// inserting a new exon into the set
	ExonSet & insert(long start, long end) {
		assert(start <= end);
		ExonPtr exonp(new Exon(start, end));
		L_(debug2) << "Trying to insert exon " << *exonp;
		VecExonp ins_exonps;
		SetExonp::iterator lower = exonps.lower_bound(exonp);
		if (exonps.size() > 0 && lower != exonps.begin()) {
			lower--;	// guarantees that (*lower)->start < exonp->start
								// or (*lower) is the beginning of the set
		}
		SetExonp::iterator itr = lower;
		while (itr != exonps.end() && (*itr)->start < exonp->end) {
			L_(debug2) << "Examining exon " << **itr;
		       	if ((*itr)->start == exonp->start) {
				// shrink the new exon
				//exonp->start = (*itr)->end;
				// shrink the old exon if the new one is shorter 
				//(*itr)->end = min((*itr)->end, exonp->end);
				//L_(debug2) << "Case 0, shrink new exon at least. New exon: " << exonp;
			        if ((*itr)->end > exonp->end) {
				        //split the old exon
				        ExonPtr exonp2(new Exon(exonp->end, (*itr)->end));
					ins_exonps.push_back(exonp2);
					(*itr)->end = exonp->end;
					exonp->start = exonp->end;
			        }
				else {
				        exonp->start = (*itr)->end;
				}
			} else if ((*itr)->start < exonp->start) {
				if ((*itr)->end > exonp->start) { // the two exons overlap
					if ((*itr)->end > exonp->end) {
						// split the old exon
						ExonPtr exonp2(new Exon(exonp->start, exonp->end));
						ExonPtr exonp3(new Exon(exonp->end, (*itr)->end));
						ins_exonps.push_back(exonp2);
						ins_exonps.push_back(exonp3);
						(*itr)->end = exonp->start;
						exonp->start = exonp->end;
						L_(debug2) << "Case 1, split old exon. New exon: " << exonp;
					} else {
						// split both exons
						ExonPtr exonp2(new Exon(exonp->start, (*itr)->end));
						ins_exonps.push_back(exonp2);
						long old_end = (*itr)->end;
						(*itr)->end = exonp->start;
						exonp->start = old_end;
						L_(debug2) << "Case 2, split both exons. New exon: " << exonp;
					}
				} else {
					// do nothing
					L_(debug2) << "Case 3, do nothing. New exon: " << exonp;
				}
			} else { // (*itr)->start > exonp->start
				if ((*itr)->end > exonp->end){
					// split both exons
					ExonPtr exonp3(new Exon(exonp->end, (*itr)->end));
					ins_exonps.push_back(exonp3);
					(*itr)->end = exonp->end;
					exonp->end = (*itr)->start;
					L_(debug2) << "Case 4, split both exons. New exon: " << exonp;
				} else {
					// split new exon
					ExonPtr exonp1(new Exon(exonp->start, (*itr)->start));
					ins_exonps.push_back(exonp1);
					exonp->start = (*itr)->end;
					L_(debug2) << "Case 5, split new exon. New exon: " << exonp;
				}
			}
			itr++;
		}
		if (exonp->start < exonp->end) {
			exonps.insert(exonp);
		}
		BOOST_FOREACH (ExonPtr const & ins_exonp, ins_exonps) {
			if (ins_exonp->start < ins_exonp->end) {
				exonps.insert(ins_exonp);
			}
		}
		L_(debug2) << "Exon set: " << exonps;
		return *this;
	}
};

class Isoforms {
public:
	string chrom;	///< The chromosome the gene is on
	string strand;	///< The strand info, "+"/"-"
	unsigned long num_exons;	///< Number of exons, N

	/// Number of all exons including ``START" and ``END", N+2
	unsigned long num_general_exons;
	unsigned long start_virtual_exon_idx;	///< The index of the ``START" exon, N
	unsigned long end_virtual_exon_idx;	///< The index of the ``END" exon, N+1

	/// Sizes of the exons, including ``START" and ``END" (both with length 0 and located at the end of the vector), a vector of size N+2
	vector<unsigned long> exon_lengths;

	/// Number of known isoforms of the gene, K1
	unsigned long num_known_isoforms;

	/// Number of all possible isoforms according to the splicing graph constructed, K2
	unsigned long num_possible_isoforms;

	/// A vector of size K1, recording the names of the known isoforms
	vector<string> known_iso_names;

	/// Probabibilties of known isoforms, an array of size K1
	double * known_iso_probs;

	/// Probabibilties of all possible isoforms, a vector of size K2
	vector<double> possible_iso_probs;

	/// K1 vectors, each recording the exon indices of a known isoform
	vector<vector<unsigned long> > known_iso_exon_indices;

	/// K1 vectors, each recording the partial total lengths of a known isoform up to (inclusive) one of its exons
	vector<vector<unsigned long> > known_iso_exon_total_lengths;

	/// K2 vectors, each recording the exon indices (excluding ``START" and ``END") of a possible isoform
	vector<vector<unsigned long> > possible_iso_exon_indices;

	/// K2 vectors, each recording the partial total lengths of a possible isoform up to (inclusive) one of its exons
	vector<vector<unsigned long> > possible_iso_exon_total_lengths;

	/// A vector of size K1, each records the index of the corresponding entry of a known isoform within the possible isoforms
	vector<unsigned long> known_iso_indices;

	Isoforms () {
		known_iso_probs = NULL;
	}

	~Isoforms () {
		delete [] known_iso_probs;
	}

	void rebuild_possible_iso_probs() {
		possible_iso_probs.resize(num_possible_isoforms, 0);
		for (unsigned long i = 0; i < num_known_isoforms; ++i) {
			possible_iso_probs[known_iso_indices[i]] = known_iso_probs[i];
		}
	}
	
	//! Constructs the information about exons and known isoforms

	/// based on the given exon set and known isoforms, also prepares the
	/// adjacency matrix (the splicing graph)
	void build(ExonSet const & exons,
			vec_gap const & iso_gaps) {
		num_exons = exons.exonps.size();
		num_known_isoforms = iso_gaps.size();
		if (num_known_isoforms > 0) {
			chrom = iso_gaps[0]->chrom;
			strand = iso_gaps[0]->strand;
		}
		known_iso_probs = new double[num_known_isoforms];
		for (unsigned long i = 0; i < num_known_isoforms; ++i) {
			known_iso_probs[i] = (double)1.0 / (double)num_known_isoforms;
		}
		start_virtual_exon_idx = num_exons;
		end_virtual_exon_idx = num_exons + 1;
		num_general_exons = num_exons + 2;
		build_isoform_array(exons, iso_gaps);
		build_splicing_graph(exons, iso_gaps);
	}

	void write_dot(ostream & os) {
		os << "digraph splicing_graph {" << endl;
		os << "\trankdir=LR;" << endl;
		for (unsigned long i = 0; i < num_general_exons; ++i) {
			if (i < start_virtual_exon_idx) {
				os << "\tnode" << i << " [label = \"Exon " << i << "\\n"
					<< exon_lengths[i] << " bp\", shape = rectangle];" << endl;
			} else if (i == start_virtual_exon_idx) {
				os << "\tnode" << i << " [label = \"START\", shape = doublecircle];" << endl;
			} else {
				os << "\tnode" << i << " [label = \"END\", shape = doublecircle];" << endl;
			}
			for (unsigned long j = 0; j < num_general_exons; ++j) {
				if (adjacency_matrix(i,j) == 1) {
					os << "\tnode" << i << " -> " << "node" << j << ";" << endl;
				}
			}
		}
		os << "}" << endl;
	}

	void enumerate_possible_isoforms() {
		vector<unsigned long> indices, total_lengths;
		for (unsigned long i = 0; i < num_exons; ++i) {
			if (adjacency_matrix(start_virtual_exon_idx, i) == 1) {
				enumerate_iso(i, indices, total_lengths);
			}
		}
		num_possible_isoforms = possible_iso_exon_indices.size();
		possible_iso_array = ublas::zero_matrix<int>(
				possible_iso_exon_indices.size(), num_exons);
		for (unsigned long i = 0;
				i < possible_iso_array.size1(); ++i) {
			BOOST_FOREACH (unsigned long idx, possible_iso_exon_indices[i]) {
				possible_iso_array(i, idx) = 1;
			}
		}
		known_iso_indices.resize(num_known_isoforms, 0);
		for (unsigned long i = 0; i < possible_iso_array.size1(); ++i) {
			for (unsigned long j = 0; j < known_iso_array.size1(); ++j) {
				bool is_same_isoform = true;
				for (unsigned long k = 0; k < num_exons; ++k){
					if (possible_iso_array(i,k) != known_iso_array(j,k)) {
						is_same_isoform = false;
						break;
					}
				}
				if (is_same_isoform) {
					known_iso_indices[j] = i;
				}
			}
		}
	}

protected:
	/// 0-1 adjacency matrix of all the exons (including ``START" and ``END", a matrix of size (N+2)*(N+2), whose (i,j) element records whether exon i connects to exon j (directional)
	ublas::matrix<int> adjacency_matrix;

	/// 0-1 matrix of size K1*N, whose (i,j) element records whether the known isoform i contains exon j
	ublas::matrix<int> known_iso_array;

	/// 0-1 matrix of size K2*N, whose (i,j) element records whether the possible isoform i contains exon j
	ublas::matrix<int> possible_iso_array;

	void build_isoform_array(ExonSet const & exons,
			vec_gap const & iso_gaps) {
		known_iso_array = ublas::zero_matrix<int>(num_known_isoforms, num_exons);
		unsigned long iso_idx = 0;
		BOOST_FOREACH (ga_ptr const & gap, iso_gaps) {
			known_iso_names.push_back(gap->gname);
			unsigned long exon_idx = 0;
			unsigned long iso_exon_idx = 0;
			BOOST_FOREACH (ExonPtr const & exonp, exons.exonps) {
				for (unsigned long i = iso_exon_idx;
						i < gap->exonCount; ++i) {
					if (exonp->is_inside(gap->exonStarts[i], gap->exonEnds[i])) {
						known_iso_array(iso_idx, exon_idx) = 1;
						iso_exon_idx = i;
						break;
					}
				}
				exon_idx++;
			}
			iso_idx++;
		}

		exon_lengths.clear();
		exon_lengths.resize(num_general_exons, 0);
		unsigned long exon_idx = 0;
		BOOST_FOREACH (ExonPtr const & exonp, exons.exonps) {
			exon_lengths[exon_idx] = exonp->end - exonp->start;
			exon_idx++;
		}

		for (unsigned long i = 0; i < num_known_isoforms; ++i) {
			vector<unsigned long> exon_indices, total_lengths;
			unsigned long total_length = 0;
			for (unsigned long j = 0; j < num_exons; ++j) {
				if (known_iso_array(i, j) == 1) {
					exon_indices.push_back(j);
					total_length += exon_lengths[j];
					total_lengths.push_back(total_length);
				}
			}
			known_iso_exon_indices.push_back(exon_indices);
			known_iso_exon_total_lengths.push_back(total_lengths);
		}
	}

	void build_splicing_graph(ExonSet const & exons,
			vec_gap const & iso_gaps) {
		// initialize adjacency matrix
		adjacency_matrix = ublas::zero_matrix<int>(
				num_general_exons, num_general_exons);
		// traverse through the isoform-exon array
		for (unsigned long i = 0; i < num_known_isoforms; ++i) {
			bool found_start = false;
			unsigned long last_exon_idx = 0;
			for (unsigned long j = 0; j < num_exons; ++j) {
				if (known_iso_array(i,j) == 1) {
					if (found_start) {
						adjacency_matrix(last_exon_idx, j) = 1;
					}
					if (!found_start) {
						adjacency_matrix(start_virtual_exon_idx,
								j) = 1;
						found_start = true;
					}
					last_exon_idx = j;
				}
			}
			adjacency_matrix(last_exon_idx,
					end_virtual_exon_idx) = 1;
		}
	}

	void enumerate_iso(unsigned long exon_idx,
			vector<unsigned long> & indices,
			vector<unsigned long> & total_lengths) {
		if (exon_idx == end_virtual_exon_idx) {
			possible_iso_exon_indices.push_back(indices);
			possible_iso_exon_total_lengths.push_back(total_lengths);
		} else {
			indices.push_back(exon_idx);
			unsigned long last_total_length = 0;
			if (total_lengths.size() > 0) {
				last_total_length = total_lengths[total_lengths.size() - 1];
			}
			total_lengths.push_back(last_total_length + exon_lengths[exon_idx]);
			for (unsigned long i = exon_idx + 1; i < num_general_exons; ++i) {
				if (adjacency_matrix(exon_idx, i) == 1) {
					enumerate_iso(i, indices, total_lengths);
				}
			}
			indices.pop_back();
			total_lengths.pop_back();
		}
	}

};	// end of class isoforms

typedef shared_ptr<Isoforms> IsoformsPtr;

#endif

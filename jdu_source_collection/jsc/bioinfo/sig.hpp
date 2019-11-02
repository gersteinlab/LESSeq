#ifndef _jsc_bioinfo_sig_hpp_included_
#define _jsc_bioinfo_sig_hpp_included_

#include <boost/config.hpp>

#include <math.h>

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include "jsc/util/interval_list.hpp"

using namespace std;
using namespace boost;
using namespace boost::lambda;

using namespace jsc::util;

namespace jsc
{

namespace bioinfo
{

//! Defines the format of a simple sig entry, without strand information.
class sig_entry
{
public:
	string chrName;
	long start;
	long end;
	double sig;
};

//! A shared_ptr to sig_entry.
/*!
 * \sa str2sigp().
 */
typedef shared_ptr<sig_entry> sig_ptr;
//! A vector of jsc::bioinfo::sig_ptr.
/*!
 * \sa load_sigps_noheader().
 */
typedef vector<sig_ptr> vec_sigp;
//! A map of (chrName, jsc::util::interval_list).
/*!
 * \sa get_map_chr_intervals().
 */
typedef map<string, interval_list<long> > map_chr_intervals;

double generalized_log2(double x) {
	if (x > 0) {
		return log2(x + 1);
	} else {
		return -log2(-x+1);
	}
}

//! Converts a string to jsc::bioinfo::sig_ptr
sig_ptr str2sigp(string const & line, bool apply_generalized_log2 = false) {
	sig_ptr sigp(new sig_entry());
	istringstream iss(line);
	iss >> sigp->chrName
		>> sigp->start
		>> sigp->end
		>> sigp->sig;
	if (apply_generalized_log2) {
		sigp->sig = generalized_log2(sigp->sig);
	}

	return sigp;
}

//! Loads a vector of jsc::bioinfo::sig_ptr from an istream (w/o the header).
vec_sigp load_sigps_noheader(istream & is)
{
	string line;

	/* read content */
	vec_sigp sigps;
	while (getline(is, line) && !is.eof())
	{
		sigps.push_back(str2sigp(line));
	}

	return sigps;
}

//! Loads a vector of jsc::bioinfo::sig_ptr from an istream (w/o the header) for a specific chrName.
void load_sigps_noheader(istream & is, vec_sigp & sigps)
{
	string line;

	/* read content */
	while (getline(is, line) && !is.eof())
	{
		sig_ptr sigp = str2sigp(line);
		sigps.push_back(str2sigp(line));
	}
}

//! Loads a vector of jsc::bioinfo::sig_ptr from an istream (w/o the header) for a specific chrName.
vec_sigp load_sigps_noheader(istream & is, string const & chrName)
{
	string line;

	/* read content */
	vec_sigp sigps;
	while (getline(is, line) && !is.eof())
	{
		sig_ptr sigp = str2sigp(line);
		if (sigp->chrName == chrName)
		{
			sigps.push_back(str2sigp(line));
		}
	}

	return sigps;
}

//! Loads a map from chrName to vector of jsc::bioinfo::sig_ptr from an istream (w/o the header)
void load_chr_sigps_map_noheader(istream & is, map<string, vec_sigp> & chr_sigps_map, bool apply_generalized_log2 = false)
{
	string line;

	/* read content */
	while (getline(is, line) && !is.eof())
	{
		sig_ptr sigp = str2sigp(line, apply_generalized_log2);
		chr_sigps_map[sigp->chrName].push_back(sigp);
	}
}

void sort_sigps(vec_sigp & sigps) {
	sort(sigps.begin(), sigps.end(),
			bind(&sig_entry::start, cref(*_1)) < bind(&sig_entry::start, cref(*_2)));
}

//! Returns all intervals whose signal values are greater/no greater than the given threshold.
/*!
 * \param sigps A vector of jsc::bioinfo::sig_ptr.
 * \param threshold The threshold.
 * \param greater_than_threshold A boolean variable specifying whether the criteria is > or <= threshold.
 * \return map_chr_intervals: a map of (chrName, jsc::util::interval_list).
 */
map_chr_intervals get_map_chr_intervals(vec_sigp const & sigps, double threshold, bool greater_than_threshold)
{
	map_chr_intervals cil_map;
	unsigned long i;
	for (i = 0; i < sigps.size(); ++i)
	{
		sig_ptr sigp = sigps[i];
		if ((greater_than_threshold && sigp->sig > threshold) ||
				(!greater_than_threshold && sigp->sig <= threshold))
		{
			cil_map[sigp->chrName].add_interval(sigp->start, sigp->end);
		}
	}

	return cil_map;
}

void get_map_chr_intervals(vec_sigp const & sigps, double threshold, bool greater_than_threshold, map_chr_intervals & cil_map)
{
	unsigned long i;
	for (i = 0; i < sigps.size(); ++i)
	{
		sig_ptr sigp = sigps[i];
		if ((greater_than_threshold && sigp->sig > threshold) ||
				(!greater_than_threshold && sigp->sig <= threshold))
		{
			cil_map[sigp->chrName].add_interval(sigp->start, sigp->end);
		}
	}
}

} /* end of bioinfo */

} /* end of jsc */

#endif

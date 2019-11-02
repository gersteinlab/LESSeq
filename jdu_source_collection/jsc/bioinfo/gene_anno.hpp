#ifndef _jsc_bioinfo_gene_anno_hpp_included_
#define _jsc_bioinfo_gene_anno_hpp_included_

#include <boost/config.hpp>

#include <math.h>

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>

#include "jsc/util/interval_list.hpp"

using namespace std;
using namespace boost;
using namespace boost::lambda;

using namespace jsc::util;

namespace jsc
{

namespace bioinfo
{

//! The format of a gene annotation entry.
class ga_entry
{
public:
	string gname;				/* 1 */
	string chrom;				/* 2 */
	string strand;				/* 3 */
	long txStart;				/* 4 */
	long txEnd;					/* 5 */
	long cdsStart;				/* 6 */
	long cdsEnd;				/* 7 */
	unsigned long exonCount;	/* 8 */
	vector<long> exonStarts;	/* 9 */
	vector<long> exonEnds;		/* 10 */
	// added to get the alternative name of the gene
	double score;
	string name2;
};

typedef shared_ptr<ga_entry> ga_ptr;
typedef vector<ga_ptr> vec_gap;

typedef map<string, ga_ptr> map_gname_gap;

//! Converts a ga line to ga_entry
ga_ptr str2gap(string const & line, bool additional_info = true)
{
	ga_ptr gap(new ga_entry());
	istringstream iss(line);
	string strExonStarts, strExonEnds;

	iss >> gap->gname			/* 1 */
		>> gap->chrom				/* 2 */
		>> gap->strand			/* 3 */
		>> gap->txStart			/* 4 */
		>> gap->txEnd				/* 5 */
		>> gap->cdsStart		/* 6 */
		>> gap->cdsEnd			/* 7 */
		>> gap->exonCount		/* 8 */
		>> strExonStarts		/* 9 */
		>> strExonEnds;			/* 10 */
	if (additional_info) {
		iss	>> gap->score				/* 11 */
			>> gap->name2;			/* 12 */
	}

	/* deal with exonStarts and exonEnds */
	typedef tokenizer<char_separator<char> > tok;
	char_separator<char> sep(",");
	tok tExonStarts(strExonStarts, sep);
	tok tExonEnds(strExonEnds, sep);

	for_each(tExonStarts.begin(), tExonStarts.end(),
			bind(&vector<long>::push_back,
				ref(gap->exonStarts),
				bind(atol,
					bind(&string::c_str, cref(_1)))));

	for_each(tExonEnds.begin(), tExonEnds.end(),
			bind(&vector<long>::push_back,
				ref(gap->exonEnds),
				bind(atol,
					bind(&string::c_str, cref(_1)))));

	return gap;
}

//! Loads a vector of ga_ptr from an istream (w/o the ga header)
vec_gap load_gaps_noheader(istream & is, bool additional_info = true)
{
	string line;

	/* read content */
	vec_gap gaps;
	while (getline(is, line) && !is.eof())
	{
		gaps.push_back(str2gap(line, additional_info));
	}

	return gaps;
}

} /* end of bioinfo */

} /* end of jsc */

#endif


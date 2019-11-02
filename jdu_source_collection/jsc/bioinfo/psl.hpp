#ifndef _jsc_bioinfo_psl_hpp_included_
#define _jsc_bioinfo_psl_hpp_included_

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

/*!
 * the format of a psl_entry.
 * according to the description at:
 * http://genome.ucsc.edu/FAQ/FAQformat
 */
class psl_entry
{
public:
	double identity;	// identity score
	double score;	// psl score
	long matches;		/* 1 */
	long misMatches;	/* 2 */
	long repMatches;	/* 3 */
	long nCount;		/* 4 */
	long qNumInsert;	/* 5 */
	long qBaseInsert;	/* 6 */
	long tNumInsert;	/* 7 */
	long tBaseInsert;	/* 8 */
	string strand;		/* 9 */
	string qName;		/* 10 */
	long qSize;			/* 11 */
	long qStart;		/* 12 */
	long qEnd;			/* 13 */
	string tName;		/* 14 */
	long tSize;			/* 15 */
	long tStart;		/* 16 */
	long tEnd;			/* 17 */
	unsigned long blockCount;	/* 18 */
	vector<long> blockSizes;	/* 19 */
	vector<long> qStarts;		/* 20 */
	vector<long> tStarts;		/* 21 */
};

typedef shared_ptr<psl_entry> psl_ptr;
typedef vector<psl_ptr> vec_pslp;
typedef map<psl_ptr, double> map_pslp_fitness;

typedef multimap<string, psl_ptr> map_qname_pslp;
typedef set<string> set_qname;

/*!
 * print a psl entry
 */
void print_pslp(psl_ptr const & pslp, ostream & os)
{
	os << pslp->matches		/* 1 */
		<< "\t" << pslp->misMatches		/* 2 */
		<< "\t" << pslp->repMatches		/* 3 */
		<< "\t" << pslp->nCount			/* 4 */
		<< "\t" << pslp->qNumInsert		/* 5 */
		<< "\t" << pslp->qBaseInsert	/* 6 */
		<< "\t" << pslp->tNumInsert		/* 7 */
		<< "\t" << pslp->tBaseInsert	/* 8 */
		<< "\t" << pslp->strand			/* 9 */
		<< "\t" << pslp->qName			/* 10 */
		<< "\t" << pslp->qSize			/* 11 */
		<< "\t" << pslp->qStart			/* 12 */
		<< "\t" << pslp->qEnd			/* 13 */
		<< "\t" << pslp->tName			/* 14 */
		<< "\t" << pslp->tSize			/* 15 */
		<< "\t" << pslp->tStart			/* 16 */
		<< "\t" << pslp->tEnd			/* 17 */
		<< "\t" << pslp->blockCount	/* 18 */
		<< "\t";
	// blockSizes						/* 19 */
	for (unsigned int i = 0; i < pslp->blockCount; ++i)
	{
		os << pslp->blockSizes[i] << ",";
	}
	os << "\t";
	// qStarts							/* 20 */
	for (unsigned int i = 0; i < pslp->blockCount; ++i)
	{
		os << pslp->qStarts[i] << ",";
	}
	os << "\t";
	// tStarts							/* 21 */
	for (unsigned int i = 0; i < pslp->blockCount; ++i)
	{
		os << pslp->tStarts[i] << ",";
	}
	os << endl;
}

/*!
 * convert a psl line to psl_entry
 */
psl_ptr str2pslp(string const & line)
{
	psl_ptr pslp(new psl_entry());
	istringstream iss(line);
	string strBlockSizes, strQStarts, strTStarts;

	iss >> pslp->matches		/* 1 */
		>> pslp->misMatches		/* 2 */
		>> pslp->repMatches		/* 3 */
		>> pslp->nCount			/* 4 */
		>> pslp->qNumInsert		/* 5 */
		>> pslp->qBaseInsert	/* 6 */
		>> pslp->tNumInsert		/* 7 */
		>> pslp->tBaseInsert	/* 8 */
		>> pslp->strand			/* 9 */
		>> pslp->qName			/* 10 */
		>> pslp->qSize			/* 11 */
		>> pslp->qStart			/* 12 */
		>> pslp->qEnd			/* 13 */
		>> pslp->tName			/* 14 */
		>> pslp->tSize			/* 15 */
		>> pslp->tStart			/* 16 */
		>> pslp->tEnd			/* 17 */
		>> pslp->blockCount		/* 18 */
		>> strBlockSizes		/* 19 */
		>> strQStarts			/* 20 */
		>> strTStarts;			/* 21 */

	/* deal with blockSizes, qStarts, and TStarts */
	typedef tokenizer<char_separator<char> > tok;
	char_separator<char> sep(",");
	tok tBlockSizes(strBlockSizes, sep);
	tok tQStarts(strQStarts, sep);
	tok tTStarts(strTStarts, sep);

	for_each(tBlockSizes.begin(), tBlockSizes.end(),
			bind(&vector<long>::push_back,
				ref(pslp->blockSizes),
				bind(atol,
					bind(&string::c_str, _1))));
			
	for_each(tQStarts.begin(), tQStarts.end(),
			bind(&vector<long>::push_back,
				ref(pslp->qStarts),
				bind(atol,
					bind(&string::c_str, _1))));

	for_each(tTStarts.begin(), tTStarts.end(),
			bind(&vector<long>::push_back,
				ref(pslp->tStarts),
				bind(atol,
					bind(&string::c_str, _1))));

	assert(pslp->qSize >= pslp->qEnd - pslp->qStart);
	assert(pslp->blockCount == pslp->blockSizes.size()
			&& pslp->blockCount == pslp->qStarts.size()
			&& pslp->blockCount == pslp->tStarts.size());
			
	return pslp;
}

/*!
 * load a vector of psl_ptrs from an istream (w/o the psl header)
 */
vec_pslp load_pslps_noheader(istream & is)
{
	string line;

	/* read content */
	vec_pslp pslps;
	while (getline(is, line) && !is.eof())
	{
		pslps.push_back(str2pslp(line));
	}

	return pslps;
}

/*!
 * a fitness score simulating the percent identity score defined at:
 * http://genome.ucsc.edu/FAQ/FAQblat.html
 *
 * here we assume that the psl_entry is a dna vs. dna blat result,
 * and we focus on the fitness of the query sequence.
 */
double psl_q_fitness_ucsc(psl_ptr const & pslp, bool const & consider_ins_factor = true)
{
	double badness = 0;
	double qAliSize = pslp->qEnd - pslp->qStart;
	double tAliSize = pslp->tEnd - pslp->tStart;
	double aliSize = min(qAliSize, tAliSize);
	if (aliSize <= 0) {
		return 0;
	}
	double sizeDif = qAliSize - tAliSize;
	if (sizeDif < 0) {
		sizeDif = 0;
	}
	double insertFactor = pslp->qNumInsert;
	if (!consider_ins_factor) {
		insertFactor = 0;
	}
	long total = pslp->matches + pslp->repMatches + pslp->misMatches;
	assert(total != 0);
	badness = (1000 * ((double)pslp->misMatches + insertFactor + round(3 * log(1 + sizeDif)))) / (double)total;

	return (100.0 - (long)badness * 0.1);
}

double psl_score_ucsc(psl_ptr const & pslp) {
	return pslp->matches +
		(pslp->repMatches >> 1) -
		pslp->misMatches -
		pslp->qNumInsert -
		pslp->tNumInsert;
}

double psl_identity(psl_ptr const & pslp) {
	return (double)(pslp->matches + pslp->repMatches) / (double)(pslp->matches + pslp->repMatches + pslp->misMatches);
}

/*!
 * a fitness score that is similar to the percent identity score
 * defined at:
 * http://genome.ucsc.edu/FAQ/FAQblat.html
 *
 * here we assume that the psl_entry is a dna vs. dna blat result,
 * and we focus on the fitness of the query sequence.
 */
double psl_q_fitness(psl_ptr const & pslp)
{
	double badness = 0;
	double sizeDif = abs((pslp->tEnd - pslp->tStart) - (pslp->qEnd - pslp->qStart))
		+ abs(pslp->qSize - (pslp->qEnd - pslp->qStart));
	double insertFactor = pslp->qNumInsert + pslp->tNumInsert;
	long total = pslp->matches + pslp->repMatches + pslp->misMatches;
	assert(total != 0);
	badness = (1000 * ((double)pslp->misMatches + (double)insertFactor + 3 * log(1 + sizeDif))) / (double)total;

	return (100 - badness * 0.1);
}

/*!
 * Note that the returned psl_entry may not be fully filled:
 * only strand, tName, tStart, tEnd, blockCount, blockSizes, tStarts
 * will be filled.
 */
psl_ptr combine_forward_reverse_reads(psl_ptr const & forward_pslp, psl_ptr const & reverse_pslp, bool & combined)
{
	assert(forward_pslp->strand != reverse_pslp->strand);
	psl_ptr combined_pslp(new psl_entry());

	combined_pslp->qName = forward_pslp->qName + "::" + reverse_pslp->qName;
	combined_pslp->strand = forward_pslp->strand;
	combined_pslp->tName = forward_pslp->tName;

	combined = true;
	interval_list<long> il_f, il_r;
	il_f.add_starts_sizes(forward_pslp->tStarts, forward_pslp->blockSizes);
	il_r.add_starts_sizes(reverse_pslp->tStarts, reverse_pslp->blockSizes);

	if (!il_f.coverage_overlap(il_r))
	{
		combined = false;
	}

	interval_list<long> il;
	il.add_interval_list(il_f);
	il.add_interval_list(il_r);

	unsigned long n = il.get_starts().size();
	assert(n > 0);
	combined_pslp->tStart = il.get_starts()[0];
	combined_pslp->tEnd = il.get_ends()[n - 1];
	combined_pslp->blockCount = n;
	for (unsigned long i = 0; i < n; i++)
	{
		combined_pslp->tStarts.push_back(il.get_starts()[i]);
		combined_pslp->blockSizes.push_back(il.get_ends()[i] - il.get_starts()[i]);
	}

	return combined_pslp;
}

/*!
 * Note that the returned psl_entry may not be fully filled:
 * only qName, strand, tName, tStart, tEnd, blockCount, blockSizes, tStarts
 * will be filled.
 */
psl_ptr fill_in_gaps(psl_ptr const & pslp, double const & threshold)
{
	psl_ptr processed_pslp(new psl_entry());
	processed_pslp->qName = pslp->qName;
	processed_pslp->strand = pslp->strand;
	processed_pslp->tName = pslp->tName;
	interval_list<long> il;
	il.add_starts_sizes(pslp->tStarts, pslp->blockSizes);
	il.fill_in_gaps(threshold);

	unsigned long n = il.get_starts().size();
	assert(n > 0);
	processed_pslp->tStart = il.get_starts()[0];
	processed_pslp->tEnd = il.get_ends()[n - 1];
	processed_pslp->blockCount = n;
	for (unsigned long i = 0; i < n; i++)
	{
		processed_pslp->tStarts.push_back(il.get_starts()[i]);
		processed_pslp->blockSizes.push_back(il.get_ends()[i] - il.get_starts()[i]);
	}

	return processed_pslp;
}

} /* end of bioinfo */

} /* end of jsc */

#endif

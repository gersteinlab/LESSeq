#ifndef _jsc_bioinfo_fasta_hpp_included_
#define _jsc_bioinfo_fasta_hpp_included_

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

using namespace std;
using namespace boost;
using namespace boost::lambda;

namespace jsc
{

namespace bioinfo
{

//! Defines the format of a simple fasta entry, without comments.
class fa_entry
{
public:
	//! Description of the sequence.
	string seq_desc;
	unsigned long seq_length;
	shared_ptr<char> seqp;
};

//! A shared_ptr to fa_entry
typedef shared_ptr<fa_entry> fa_ptr;
//! A vector of fa_ptrs
typedef vector<fa_ptr> vec_fap;
//! desc to fap map
typedef map<string, fa_ptr> map_desc_fap;

//! Loads a fasta_entry from an istream.
fa_ptr load_fap(istream & is)
{
	vector<string> vec_string;
	fa_ptr fap(new fa_entry());
	string line;
	char c;
	bool desc_loaded = false;
	unsigned long seq_length = 0;
	while (is.peek() != EOF)
	{
		c = is.peek();
		if (c == '>' && !desc_loaded)
		{
			// assign description
			getline(is, line);
			fap->seq_desc = line;
			desc_loaded = true;
			continue;
		}
		else if (c == '>' && desc_loaded)
		{
			break;
		}
		getline(is, line);
		seq_length += line.size();
		vec_string.push_back(line);
	}
	// assign seq_length
	fap->seq_length = seq_length;
	// allocate memory for fap->seq
	shared_ptr<char> p(new char[fap->seq_length]);
	unsigned long cur_pos = 0;
	unsigned long i;
	for (i = 0; i < vec_string.size(); ++i)
	{
		memcpy(p.get() + cur_pos, vec_string[i].c_str(), vec_string[i].size());
		cur_pos += vec_string[i].size();
	}
	fap->seqp = p;

	return fap;
}

//! Loads fa_ptrs from an istream.
vec_fap load_faps_noheader(istream & is)
{
	/* read content */
	vec_fap faps;
	while (!is.eof())
	{
		faps.push_back(load_fap(is));
	}

	return faps;
}

} /* end of bioinfo */

} /* end of jsc */

#endif

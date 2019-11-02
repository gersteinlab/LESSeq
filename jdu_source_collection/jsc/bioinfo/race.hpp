#ifndef _jsc_bioinfo_race_hpp_included_
#define _jsc_bioinfo_race_hpp_included_

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

class race_design_entry
{
public:
	string qname;
	string primer_name;
	string pcr_direction;	/* 5' or 3' */
	string read_direction;	/* forward/backward_read */
};

typedef shared_ptr<race_design_entry> race_ptr;
typedef vector<race_ptr> vec_racep;
typedef map<string, race_ptr> map_qname_racep;
typedef multimap<string, race_ptr> map_primer_racep;
typedef set<string> set_primer;

race_ptr str2racep(string const & line)
{
	race_ptr racep(new race_design_entry());
	istringstream iss(line);
	iss >> racep->qname
		>> racep->primer_name
		>> racep->pcr_direction
		>> racep->read_direction;

	assert(racep->pcr_direction == "5'" || racep->pcr_direction == "3'");
	assert(racep->read_direction == "forward_read" || racep->read_direction == "reverse_read");

	return racep;
}

vec_racep load_raceps_noheader(istream & is)
{
	string line;
	vec_racep raceps;
	while (getline(is, line) && !is.eof())
	{
		raceps.push_back(str2racep(line));
	}

	return raceps;
}

} /* end of bioinfo */

} /* end of jsc */

#endif

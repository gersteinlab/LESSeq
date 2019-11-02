#ifndef _jsc_bioinfo_bartxt_hpp_included_
#define _jsc_bioinfo_bartxt_hpp_included_

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

//! Defines the format of a bar entry, without strand information.
class bar_entry
{
public:
	long pos;
	double sig;
};

//! A shared_ptr to bar_entry.
typedef shared_ptr<bar_entry> bar_ptr;
//! A vector of jsc::bioinfo::bar_ptr.
typedef vector<bar_ptr> vec_barp;
//! A map of (chrName, vector of bar_entries).
typedef map<string, vec_barp> map_chr_barps;

//! Loads map_chr_barps from an istream (w/o the header).
map_chr_barps load_barps(istream & is)
{
	string line;
	/* read content */
	map_chr_barps chr_barps_map;
	string chrName = "";
	string prefix = "# Name";
	while (getline(is, line) && !is.eof())
	{
		if (line.length() == 0)
		{
			continue;
		}
		else if (line.compare(0, prefix.size(), prefix) == 0)
		{
			// load chrName
			istringstream iss(line);
			string s1, s2; //``#" and ``Name"
			iss >> s1
				>> s2
				>> chrName;
		}
		else if (line[0] == '#')
		{
			continue;
		}
		else
		{
			// load bar_entry
			bar_ptr barp(new bar_entry());
			istringstream iss(line);
			iss >> barp->pos
				>> barp->sig;
			chr_barps_map[chrName].push_back(barp);
		}
	}

	// sort all the vectors in the map
	map_chr_barps::iterator cbm_itr;
	for (cbm_itr = chr_barps_map.begin(); cbm_itr != chr_barps_map.end(); ++cbm_itr)
	{
		sort((*cbm_itr).second.begin(), (*cbm_itr).second.end(),
				bind(&bar_entry::pos, *_1) < bind(&bar_entry::pos, *_2));
	}
	return chr_barps_map;
}

} /* end of bioinfo */

} /* end of jsc */

#endif

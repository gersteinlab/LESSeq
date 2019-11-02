#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include "../common/read.h"

using namespace std;
using namespace boost;
using namespace boost::lambda;

int syntax_err() {
	L_(error) << "Usage:\nclassify\n	log_level(0,1,2,...) proj_name out_prefix\n	isoform_format isoforms_path g2i_format g2i_path gene_begin_idx gene_end_idx\n";
	return 1;
}

int lexical_cast_err() {
	L_(error) << "Lexical_cast error when converting arguments to numeric values";
	return 1;
}

int unknown_format_err(string const & format) {
	L_(error) << "Unknown file format error: " << format;
	return 1;
}

int unknown_readtype_err(string const & type) {
	L_(error) << "Unknown read type error: " << type;
	return 1;
}

class FileFormats {
	public:
		static const string UCSC_GENE2ISOFORM;
		static const string LH_GENE_TXT;
};
const string FileFormats::UCSC_GENE2ISOFORM = "UCSC_GENE2ISOFORM";
const string FileFormats::LH_GENE_TXT = "LH_GENE_TXT";

typedef tuple<string, string, interval_list<long> > gene_info;


int main(int argc, char * argv[])
{
	typedef tokenizer<char_separator<char> > tok;
	// initialize the logging system
	Log::ReportingLevel = LogLevel::info;
	
	// get parameters
	if (argc < 10) {
		return syntax_err();
	}
	int argi = 1;
	long log_level = lexical_cast<long>(argv[argi++]);
	Log::ReportingLevel = log_level;
	string proj_name = argv[argi++];
	string out_prefix = argv[argi++];

	// load isoform parameters
	string isoform_format = argv[argi++];
	string isoforms_path = argv[argi++];
	string g2i_format = argv[argi++];
	string g2i_path = argv[argi++];
	unsigned long gene_begin_idx, gene_end_idx;
	try {
		gene_begin_idx = lexical_cast<unsigned long>(argv[argi++]);
		gene_end_idx = lexical_cast<unsigned long>(argv[argi++]);

	} catch (bad_lexical_cast &) {
		return lexical_cast_err();
	}

	// global variables
	ifstream ifs;
	string line;

	// load isoforms
	L_(info) << "Loading isoforms...";
	ifs.clear();
	ifs.open(isoforms_path.c_str());
	assert(ifs.is_open() && !ifs.eof());
	vec_gap iso_gaps;
	if (isoform_format == FileFormats::LH_GENE_TXT) {
		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			ga_ptr gap(new ga_entry());

			string strExonStarts, strExonEnds;
			iss >> gap->gname >> gap->chrom >> gap->strand
				>> gap->txStart >> gap->txEnd
				>> gap->exonCount >> strExonStarts >> strExonEnds;
			gap->cdsStart = gap->txStart;
			gap->cdsEnd = gap->txEnd;

			/* deal with exonStarts and exonEnds */
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

			iso_gaps.push_back(gap);
		}
	} else {
		return unknown_format_err(isoform_format);
	}
	ifs.close();
	map<string, ga_ptr> iname2gap;
	BOOST_FOREACH(ga_ptr const & gap, iso_gaps) {
		iname2gap[gap->gname] = gap;
	}
	L_(info) << "Loaded " << iso_gaps.size() << " isoforms";

	L_(info) << "Loading gene->isoform mappings ...";
	ifs.clear();
	ifs.open(g2i_path.c_str());
	assert(ifs.is_open() && !ifs.eof());
	map<string, vector<ga_ptr> > gname2isogaps;
	set<string> gname_set;
	if (g2i_format == FileFormats::UCSC_GENE2ISOFORM) {
		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			string gname, iname;
			iss >> gname >> iname;
			gname_set.insert(gname);
			gname2isogaps[gname].push_back(iname2gap[iname]);
		}
	} else {
		return unknown_format_err(g2i_format);
	}
	ifs.close();
	L_(info) << "Loaded " << gname_set.size() << " genes";

	L_(info) << "Selecting genes: from #" << gene_begin_idx
		<< " to #" << gene_end_idx << ", ignoring single-isoform ones";
	vector<string> gnames;
	unsigned long count = 0;
	BOOST_FOREACH(string const & gname, gname_set) {
		if (count >= gene_begin_idx
				&& count < gene_end_idx) {
			if (gname2isogaps[gname].size() >= 2) 	// remove to report rpkm values for genes with single isoform as well
			{
				gnames.push_back(gname);
			}
		}
		count++;
	}
	L_(info) << "Selected " << gnames.size() << " gene(s)";
	L_(debug) << "Gene names:";
	BOOST_FOREACH(string const & gname, gnames) {
		L_(debug) << gname;
		BOOST_FOREACH(ga_ptr const gap, gname2isogaps[gname]) {
			L_(debug) << "  " << gap->gname;
		}
	}

	// For each gene, build its isoform structure...
	L_(info) << "Building isoform structures for the "
		<< gnames.size() << " selected gene(s)";
	map<string, interval_list<long> > gname2il;
	map<string, ExonSet> gname2exons;
	map<string, IsoformsPtr> gname2isop;
	BOOST_FOREACH(string const & gname, gnames) {
		// get isoforms
		vector<ga_ptr> const & iso_gaps = gname2isogaps[gname];
		// generate exon set
		ExonSet exons;
		IsoformsPtr isop(new Isoforms);
		BOOST_FOREACH(ga_ptr const & gap, iso_gaps) {
			for (unsigned long i = 0; i < gap->exonCount; ++i) {
				exons.insert(gap->exonStarts[i], gap->exonEnds[i]);
				gname2il[gname].add_interval(
						gap->exonStarts[i], gap->exonEnds[i]);
			}
		}
		gname2exons[gname] = exons;
		// build isoforms
		isop->build(exons, iso_gaps);
		gname2isop[gname] = isop;
		
        ofstream tmp_ofs;
        tmp_ofs.clear();
        string tmp_path = out_prefix + gname + ".matrix";
        tmp_ofs.open(tmp_path.c_str());
		assert(tmp_ofs.is_open() && !tmp_ofs.eof());

    	tmp_ofs << isop->chrom << "\t" << isop->strand << "\t" << exons.exonps << "\n";

		ublas::matrix<int> tmp_iso_array = ublas::zero_matrix<int>(isop->num_known_isoforms, isop->num_exons);
        unsigned long tmp_iso_idx = 0;
		vector<string> tmp_known_iso_names;
        BOOST_FOREACH (ga_ptr const & gap, iso_gaps) {
		  	tmp_known_iso_names.push_back(gap->gname);
		  	unsigned long tmp_exon_idx = 0;
		  	unsigned long tmp_iso_exon_idx = 0;
		  	BOOST_FOREACH (ExonPtr const & exonp, exons.exonps) {
		    	for (unsigned long i = tmp_iso_exon_idx;
			 					i < gap->exonCount; ++i) {
		      		if (exonp->is_inside(gap->exonStarts[i], gap->exonEnds[i])) {
						tmp_iso_array(tmp_iso_idx, tmp_exon_idx) = 1;
						tmp_iso_exon_idx = i;
						break;
		      		}
		    	}
                tmp_ofs << tmp_iso_array(tmp_iso_idx, tmp_exon_idx) << "\t";
		   		tmp_exon_idx++;
		 	}
        	tmp_ofs << "\n";
		 	tmp_iso_idx++;
    	}
	
		// no need to enumerate all possible isoforms for
		// the purpose of quantifying known isoforms
		// isop->enumerate_possible_isoforms();
		// isop->rebuild_possible_iso_probs();
	}
	L_(info) << "Built isoform structures for the "
		<< gnames.size() << " selected gene(s)";

	return 0;
}

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
	L_(error) << "Usage:\nsolve\n	log_level(0,1,2,...) proj_name out_prefix\n	isoform_format isoforms_path g2i_format g2i_path gene_begin_idx gene_end_idx\n  (read_format read_type expected_read_length reads_path total_read_bases)+";
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
		static const string UCSC_GENE_TXT;
		static const string UCSC_GENE2ISOFORM;
		static const string UCSC_GFF;
		static const string UCSC_BED;
		static const string LH_GENE_TXT;
		static const string WORMBASE_GFF2;
		static const string WORMBASE_GENE2ISOFORMS;
		static const string WORMBASE_GFF3;
		static const string MRF_SINGLE;
		static const string GENELETS_GFF3;
};
const string FileFormats::UCSC_GENE_TXT = "UCSC_GENE_TXT";
const string FileFormats::UCSC_GENE2ISOFORM = "UCSC_GENE2ISOFORM";
const string FileFormats::UCSC_GFF = "UCSC_GFF";
const string FileFormats::UCSC_BED = "UCSC_BED";
const string FileFormats::LH_GENE_TXT = "LH_GENE_TXT";
const string FileFormats::WORMBASE_GFF2 = "WORMBASE_GFF2";
const string FileFormats::WORMBASE_GENE2ISOFORMS= "WORMBASE_GENE2ISOFORMS";
const string FileFormats::WORMBASE_GFF3 = "WORMBASE_GFF3";
const string FileFormats::MRF_SINGLE = "MRF_SINGLE";
const string FileFormats::GENELETS_GFF3 = "GENELETS_GFF3";

typedef tuple<string, string, interval_list<long> > read_info;
typedef tuple<string, string, interval_list<long> > gene_info;

class rinfo {
	public:
		string name;
		string chrom;
		string strand;
		long start;
		long end;
};

typedef shared_ptr<rinfo> rinfo_ptr;

class RinfopComp {
	public:
		bool operator() (rinfo_ptr const & a, rinfo_ptr const & b) {
			if (a->chrom == b-> chrom) {
				if (a->start == b->start) {
					if (a->end == b->end) {
						if (a->strand == b->strand) {
							return (a->name < b->name);
						} else {
							return (a->strand < b->strand);
						}
					} else {
						return (a->end < b->end);
					}
				} else {
					return (a->start < b->start);
				}
			} else {
				return (a->chrom < b->chrom);
			}
		}
};
typedef set<rinfo_ptr, RinfopComp> set_rinfop;

int main(int argc, char * argv[])
{
	typedef tokenizer<char_separator<char> > tok;
	// initialize the logging system
	Log::ReportingLevel = LogLevel::info;
	
	// get parameters
	if (argc < 15) {
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
	vector<string> read_formats, read_types, reads_paths;
	vector<double> total_read_bases_vec;
	vector<unsigned long> expected_read_lengths;
	try {
		gene_begin_idx = lexical_cast<unsigned long>(argv[argi++]);
		gene_end_idx = lexical_cast<unsigned long>(argv[argi++]);

		// load read parameters
		while (argi < argc) {
			if (argc - argi < 5) {
				return syntax_err();
			}
			read_formats.push_back(argv[argi++]);
			read_types.push_back(argv[argi++]);
			expected_read_lengths.push_back(
					lexical_cast<unsigned long>(argv[argi++]));
			reads_paths.push_back(argv[argi++]);
			total_read_bases_vec.push_back(
					lexical_cast<double>(argv[argi++]));
		}
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
	if (isoform_format == FileFormats::UCSC_GENE_TXT) {
		iso_gaps = load_gaps_noheader(ifs, false);
	} else if (isoform_format == FileFormats::GENELETS_GFF3) {
		map<string, gene_info> iname2ginfo;
		set<string> inames;

		getline(ifs, line); // skip first line
		getline(ifs, line); // skip second line
		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			string type;
			long start, end;
			string infos;
			string chr, tmp, strand;
			iss >> chr >> tmp >> type;
			if (type == "exon") {
				iss >> start >> end
					>> tmp >> strand >> tmp >> infos;
				char_separator<char> sep(";");
				tok tInfo(infos, sep);
				BOOST_FOREACH (string const & info, tInfo) {
					if (info.size() > 7 &&
							info.substr(0, 7) == "Parent=") {
						string local_inames = info.substr(7);
						char_separator<char> sep2(",");
						tok tIname(local_inames, sep2);
						BOOST_FOREACH (string const & iname, tIname) {
							inames.insert(iname);
							iname2ginfo[iname].get<0>() = "chr" + chr;
							iname2ginfo[iname].get<1>() = strand;
							iname2ginfo[iname].get<2>().add_interval(start - 1, end);
						}
					}
				}
			}
		}
		BOOST_FOREACH(string const & iname, inames) {
			ga_ptr gap(new ga_entry());
			gap->gname = iname;
			gap->chrom = iname2ginfo[iname].get<0>();
			gap->strand = iname2ginfo[iname].get<1>();
			interval_list<long> const & il = iname2ginfo[iname].get<2>();
			gap->exonCount = il.get_num_intervals();
			gap->exonStarts = il.get_starts();
			gap->exonEnds = il.get_ends();
			iso_gaps.push_back(gap);
		}
	} else if (isoform_format == FileFormats::LH_GENE_TXT) {
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
	} else if (isoform_format == FileFormats::UCSC_GFF) {
		map<string, gene_info> iname2ginfo;
		set<string> inames;

		getline(ifs, line); // skip first line
		getline(ifs, line); // skip second line

		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			long start, end;
			string iname, chr, tmp, strand;
			iss >> chr >> tmp >> tmp >> start >> end
				>> tmp >> strand >> tmp >> iname;
			chr = chr;
			trim_if(iname, is_any_of("\""));
			inames.insert(iname);
			iname2ginfo[iname].get<0>() = chr;
			iname2ginfo[iname].get<1>() = strand;
			iname2ginfo[iname].get<2>().add_interval(start - 1, end);
		}
		BOOST_FOREACH(string const & iname, inames) {
			ga_ptr gap(new ga_entry());
			gap->gname = iname;
			gap->chrom = iname2ginfo[iname].get<0>();
			gap->strand = iname2ginfo[iname].get<1>();
			interval_list<long> const & il = iname2ginfo[iname].get<2>();
			gap->exonCount = il.get_num_intervals();
			gap->exonStarts = il.get_starts();
			gap->exonEnds = il.get_ends();
			iso_gaps.push_back(gap);
		}
	} else if (isoform_format == FileFormats::WORMBASE_GFF2) {
		map<string, gene_info> iname2ginfo;
		set<string> inames;
		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			long start, end;
			string iname, chr, tmp, strand;
			iss >> chr >> tmp >> tmp >> start >> end
				>> tmp >> strand >> tmp >> tmp >> iname;
			chr = "chr" + chr;
			trim_if(iname, is_any_of("\""));
			inames.insert(iname);
			iname2ginfo[iname].get<0>() = chr;
			iname2ginfo[iname].get<1>() = strand;
			iname2ginfo[iname].get<2>().add_interval(start - 1, end);
		}
		BOOST_FOREACH(string const & iname, inames) {
			ga_ptr gap(new ga_entry());
			gap->gname = iname;
			gap->chrom = iname2ginfo[iname].get<0>();
			gap->strand = iname2ginfo[iname].get<1>();
			interval_list<long> const & il = iname2ginfo[iname].get<2>();
			gap->exonCount = il.get_num_intervals();
			gap->exonStarts = il.get_starts();
			gap->exonEnds = il.get_ends();
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
	} else if (g2i_format == FileFormats::WORMBASE_GENE2ISOFORMS) {
		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			string gname, inames;
			iss >> gname >> inames;
			gname_set.insert(gname);
			char_separator<char> sep(";");
			tokenizer<char_separator<char> > tokens(inames, sep);
			BOOST_FOREACH(string const & iname, tokens) {
				gname2isogaps[gname].push_back(iname2gap[iname]);
			}
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
			//if (gname2isogaps[gname].size() >= 2) {	// removed to report rpkm values for genes with single isoform as well
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
	// also, keep track of the genomic regions covered by these genes
	// so that the reads can be filtered when they are loaded at a
	// later stage
	L_(info) << "Building isoform structures for the "
		<< gnames.size() << " selected gene(s)";
	map<string, interval_list<long> > covered_regions;
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
				covered_regions[gap->chrom].add_interval(
						gap->exonStarts[i], gap->exonEnds[i]);
				gname2il[gname].add_interval(
						gap->exonStarts[i], gap->exonEnds[i]);
			}
		}
		gname2exons[gname] = exons;
		// build isoforms
		isop->build(exons, iso_gaps);
		// no need to enumerate all possible isoforms for
		// the purpose of quantifying known isoforms
		// isop->enumerate_possible_isoforms();
		// isop->rebuild_possible_iso_probs();
		gname2isop[gname] = isop;
	}
	L_(info) << "Built isoform structures for the "
		<< gnames.size() << " selected gene(s)";

	unsigned long M = reads_paths.size();
	L_(info) << "Loading the reads from " << M
		<< " sampling method(s)";
	vector<set<string> > rnames_vec;
	vector<map<string, read_info> > rname2rinfo_vec;
	for (unsigned long m = 0; m < M; ++m) {
		L_(info) << "Loading the reads from sampling method #" << m;
		string read_format = read_formats[m];
		string reads_path = reads_paths[m];
		// load the reads
		set<string> rnames;
		map<string, read_info> rname2rinfo;
		rnames_vec.push_back(rnames);
		rname2rinfo_vec.push_back(rname2rinfo);
		ifs.clear();
		ifs.open(reads_path.c_str());
		assert(ifs.is_open() && !ifs.eof());
		if (read_format == FileFormats::UCSC_GFF) {
			getline(ifs, line); // skip first line
			getline(ifs, line); // skip second line
			while (getline(ifs, line) && !ifs.eof()) {
				istringstream iss(line);
				long start, end;
				string rname, chr, tmp, strand;
				iss >> chr >> tmp >> tmp >> start >> end
					>> tmp >> strand >> tmp >> rname;
				if (covered_regions[chr].contains_interval(start - 1, end)) {
					rnames_vec[m].insert(rname);
					rname2rinfo_vec[m][rname].get<0>() = chr;
					rname2rinfo_vec[m][rname].get<1>() = strand;
					rname2rinfo_vec[m][rname].get<2>().add_interval(start - 1, end);
				}
			}
		} else if (read_format == FileFormats::MRF_SINGLE) {
			string::size_type old_colon, colon;
			string chr, strand, rname;
			long start, end;
			getline(ifs, line); // skip first line
			unsigned long line_num = 0;
			while (getline(ifs, line) && !ifs.eof()) {
				line_num++;

				if ("#" == line.substr(0, 1) || "AlignmentBlocks" == line) {
					continue;
				}

				try {
					std::stringstream ss;
					ss << "read-" << line_num;
					rname = ss.str();

					string::size_type last_comma = 0;

					while (last_comma != string::npos) {
						colon = line.find(':', last_comma);	// end of col #1
						// col #1: chr
						chr = line.substr(
								last_comma == 0 ? 0 : last_comma + 1,
								colon - (last_comma == 0 ? 0 : last_comma + 1));
						old_colon = colon;
						colon = line.find(':', colon + 1);	// end of col #2
						// col #2: strand
						strand = line.substr(old_colon + 1, colon - old_colon - 1);
						old_colon = colon;
						colon = line.find(':', colon + 1);	// end of col #3
						// col #3: start
						start = lexical_cast<long>(line.substr(old_colon + 1,
									colon - old_colon - 1));
						old_colon = colon;
						colon = line.find(':', colon + 1);	// end of col #4
						// col #4: end
						end = lexical_cast<long>(line.substr(old_colon + 1,
									colon - old_colon - 1));
						if (covered_regions[chr].contains_interval(start - 1, end)) {
							rnames_vec[m].insert(rname);
							rname2rinfo_vec[m][rname].get<0>() = chr;
							rname2rinfo_vec[m][rname].get<1>() = strand;
							rname2rinfo_vec[m][rname].get<2>().add_interval(start - 1, end);
						}

						last_comma = line.find(',', colon);
					}
				} catch (bad_lexical_cast &) {
					L_(error) << "#" << line_num << ":" << line;
					return lexical_cast_err();
				}

				if (line_num % 1000000 == 0) {
					L_(info) << "Loaded " << line_num << " reads...";
				}
			}
		} else if (read_format == FileFormats::UCSC_BED) {
			long start, end, num_blocks, istart, isize;
			string chr, strand, rname, starts, sizes;
			string::size_type tabstop, old_tabstop,
				old_comma1, old_comma2, comma1, comma2;
			getline(ifs, line); // skip first line
			while (getline(ifs, line) && !ifs.eof()) {
				tabstop = line.find_first_of('\t');	// end of col #1
				// col #1: chr
				chr = line.substr(0, tabstop);
				old_tabstop = tabstop;
				tabstop = line.find('\t', tabstop + 1);	// end of col #2
				// col #2: start
				start = lexical_cast<long>(line.substr(old_tabstop + 1,
							tabstop - old_tabstop - 1));
				old_tabstop = tabstop;
				tabstop = line.find('\t', tabstop + 1);	// end of col #3
				end = lexical_cast<long>(line.substr(old_tabstop + 1,
							tabstop - old_tabstop - 1));
				if (!(covered_regions[chr].contains_interval(start, end))) {
					continue;
				}
				old_tabstop = tabstop;
				tabstop = line.find('\t', tabstop + 1);	// end of col #4
				// col #4: rname
				rname = line.substr(old_tabstop + 1, tabstop - old_tabstop - 1);
				tabstop = line.find('\t', tabstop + 1);	// end of col #5
				old_tabstop = tabstop;
				tabstop = line.find('\t', tabstop + 1);	// end of col #6
				// col #6: strand
				strand = line.substr(old_tabstop + 1, tabstop - old_tabstop - 1);
				tabstop = line.find('\t', tabstop + 1);	// end of col #7
				tabstop = line.find('\t', tabstop + 1);	// end of col #8
				tabstop = line.find('\t', tabstop + 1);	// end of col #9
				old_tabstop = tabstop;
				tabstop = line.find('\t', tabstop + 1);	// end of col #10
				// col #10: blockCount
				num_blocks = lexical_cast<long>(line.substr(old_tabstop + 1,
							tabstop - old_tabstop - 1));
				old_tabstop = tabstop;
				tabstop = line.find('\t', tabstop + 1);	// end of col #11
				// col #11: blockSizes
				sizes = line.substr(old_tabstop + 1, tabstop - old_tabstop - 1);
				old_tabstop = tabstop;
				tabstop = line.find('\t', tabstop + 1);	// end of col #12
				// col #12: blockStarts
				starts = line.substr(old_tabstop + 1, tabstop - old_tabstop - 1);
				rnames_vec[m].insert(rname);
				rname2rinfo_vec[m][rname].get<0>() = chr;
				rname2rinfo_vec[m][rname].get<1>() = strand;

				old_comma1 = old_comma2 = -1;
				for (long i = 0; i < num_blocks; ++i) {
					comma1 = starts.find(',', old_comma1 + 1);
					comma2 = sizes.find(',', old_comma2 + 1);
					istart = lexical_cast<long>(starts.substr(old_comma1 + 1,
								comma1 - old_comma1 - 1));
					isize = lexical_cast<long>(sizes.substr(old_comma2 + 1,
								comma2 - old_comma2 - 1));
					rname2rinfo_vec[m][rname].get<2>().
						add_interval(start + istart, start + istart + isize);
					old_comma1 = comma1;
					old_comma2 = comma2;
				}
			}
		} else if (read_format == FileFormats::WORMBASE_GFF3) {
			long start, end;
			string chr, strand, rinfo, rname, attrib, value;
			string::size_type tabstop, old_tabstop, old_found, found;
			unsigned long line_num = 0;
			while (getline(ifs, line) && !ifs.eof()) {
				line_num++;
				tabstop = line.find_first_of('\t');
				chr = "chr" + line.substr(0, tabstop);
				tabstop = line.find('\t', tabstop + 1);
				tabstop = line.find('\t', tabstop + 1);
				old_tabstop = tabstop;
				tabstop = line.find('\t', old_tabstop + 1);
				start = lexical_cast<long>(line.substr(old_tabstop + 1,
							tabstop - old_tabstop - 1));
				old_tabstop = tabstop;
				tabstop = line.find('\t', old_tabstop + 1);
				end = lexical_cast<long>(line.substr(old_tabstop + 1,
							tabstop - old_tabstop - 1));
				tabstop = line.find('\t', tabstop + 1);
				old_tabstop = tabstop;
				tabstop = line.find('\t', old_tabstop + 1);
				strand = line.substr(old_tabstop + 1,
						tabstop - old_tabstop - 1);
				tabstop = line.find('\t', tabstop + 1);
				old_tabstop = tabstop;
				tabstop = line.find('\t', old_tabstop + 1);
				rinfo = line.substr(old_tabstop + 1,
						tabstop - old_tabstop - 1);

				bool found_target = false;
				bool found_parent = false;
				unsigned long start2 = 0, end2 = 0;
				old_tabstop = -1;
				tabstop = rinfo.find_first_of(';');
				while (tabstop != string::npos) {
					attrib = rinfo.substr(old_tabstop + 1,
							tabstop - old_tabstop - 1);
					if (attrib.substr(0, 7) == "Target=") {
						value = attrib.substr(7);
						found = value.find_first_of(' ');
						rname = value.substr(0, found);
						found_target = true;
					} else if (attrib.substr(0, 7) == "Parent=") {
						value = attrib.substr(7);
						if (value.substr(0, 7) == "intron_") {
							found_parent = true;
							found = value.find('_', 7);
							old_found = found;
							found = value.find('_', found + 1);
							start2 = lexical_cast<unsigned long>(value.substr(
										old_found + 1, found - old_found - 1));
							old_found = found;
							found = value.find('_', found + 1);
							end2 = lexical_cast<unsigned long>(value.substr(
										old_found + 1, found - old_found - 1));
						}
					}

					old_tabstop = tabstop;
					tabstop = rinfo.find(';', old_tabstop + 1);
				}
				if (!(found_target)) {
					L_(error) << "No target description exists for read: "
						<< line;
				}
				
				if (covered_regions[chr].contains_interval(start - 1, end)) {
					rnames_vec[m].insert(rname);
					rname2rinfo_vec[m][rname].get<0>() = chr;
					rname2rinfo_vec[m][rname].get<1>() = strand;
					if (!found_parent) {
						rname2rinfo_vec[m][rname].get<2>().add_interval(start - 1, end);
					} else {
						rname2rinfo_vec[m][rname].get<2>().add_interval(start - 1, start2 - 1);
						rname2rinfo_vec[m][rname].get<2>().add_interval(end2, end);
					}
				}

				if (line_num % 1000000 == 0) {
					L_(info) << "Loaded " << line_num << " reads...";
				}
			}
		} else {
			return unknown_format_err(read_format);
		}
		ifs.close();
		L_(info) << "Sampling method #" << m
		 << ": loaded "	<< rnames_vec[m].size()
		 << " reads associated with the selected gene regions, with "
		 << total_read_bases_vec[m] << " mapped read bases in total";
	}

	// index the reads
	L_(info) << "Indexing the reads...";
	vector<set_rinfop> rinfops_vec;
	for (unsigned long m = 0; m < M; ++m) {
		set_rinfop rinfops;
		set<string> const & rnames = rnames_vec[m];
		BOOST_FOREACH(string const & rname, rnames) {
			rinfo_ptr rinfop(new rinfo());
			rinfop->name = rname;
			rinfop->chrom = rname2rinfo_vec[m][rname].get<0>();
			rinfop->strand = rname2rinfo_vec[m][rname].get<1>();
			interval_list<long> const & read_il =
				rname2rinfo_vec[m][rname].get<2>();
			rinfop->start = read_il.get_starts()[0];
			rinfop->end = read_il.get_ends()[read_il.get_num_intervals() - 1];
			rinfops.insert(rinfop);
		}
		rinfops_vec.push_back(rinfops);
	}

	// quantifying isoforms for each gene
	L_(info) << "Processing reads info for genes";
	unsigned long num_processed_genes = 0;
	BOOST_FOREACH(string const & gname, gnames) {
		L_(debug) << "Processing reads info for gene " << gname << "...";
		vector<ublas::matrix<double> > m_delta_Gs;
		IsoformsPtr isop = gname2isop[gname];
		interval_list<long> gene_il = gname2il[gname];
		long gene_start = gene_il.get_starts()[0];
		long gene_end = gene_il.get_ends()[gene_il.get_num_intervals() - 1];
		shared_ptr<Read> readp;
		vector<double> supports;
		vector<double> support_bases;
		// record the reads
		//ofstream ofs;
		//ofs.clear();
		//string gene_reads_path = out_prefix + gname + "-reads.gff";
		//ofs.open(gene_reads_path.c_str());
		//assert(ofs.is_open() && !ofs.eof());
		//ofs << "browser position " << isop->chrom
		//<< ":" << gene_start << "-" << gene_end << endl;
		for (unsigned long m = 0; m < M; ++m) {
		  //	ofs << "track name=\"" << proj_name
		  //		<< ": reads in " << gname
		  //		<< " from sampling method #" << m
		  //		<< "\" visibility=4" << endl;
			// assign read type
			string read_type = read_types[m];
			if (read_type == ReadType::TYPES[ReadType::MEDIUM_READ]) {
				shared_ptr<AccessibleReadStarts> ars(
						new AccessibleReadStarts(
							isop->known_iso_exon_indices,
							isop->known_iso_exon_total_lengths,
							isop->exon_lengths,
							expected_read_lengths[m]));
				ars->construct_ARS();
				shared_ptr<Read> readp_medium_single(
						new Read_medium_single(ars, NULL, expected_read_lengths[m]));
				readp = readp_medium_single;
			} else if (read_type == ReadType::TYPES[ReadType::SHORT_READ]) {
				shared_ptr<AccessibleReadStarts> ars(
						new AccessibleShortReadStarts(
							isop->known_iso_exon_indices,
							isop->known_iso_exon_total_lengths,
							isop->exon_lengths,
							expected_read_lengths[m],
							0));
				ars->construct_ARS();
				shared_ptr<Read> readp_short_single(
						new Read_short_single(ars, NULL, expected_read_lengths[m]));
				readp = readp_short_single;
			} else {
				return unknown_readtype_err(read_type);
			}
			// find compatiable reads
			vector<string> valid_rnames;
			set_rinfop & rinfops = rinfops_vec[m];
			rinfo_ptr gp(new rinfo());
			gp->name = gname;
			gp->chrom = isop->chrom;
			gp->strand = isop->strand;
			gp->start = gene_start;
			gp->end = gene_end;
			set_rinfop::const_iterator itr =
				rinfops.lower_bound(gp);
			double num_valid_read_bases = 0;
			while (itr != rinfops.end() && (*itr)->chrom == gp->chrom && (*itr)->start <= gp->end ) {
			  if ( (*itr)->start >= gp->start ) {
					string rname = (*itr)->name;
					interval_list<long> const & read_il =
						rname2rinfo_vec[m][rname].get<2>();
					readp->build(read_il, gname2exons[gname].exonps);
					bool compatible = false;
					for (unsigned long j = 0;
							j < isop->num_known_isoforms; ++j) {
						if (readp->is_compatible_with_iso(j)) {
							if ((double)readp->get_read_length()
									/ read_il.compute_total_length() > 0.98) {
								valid_rnames.push_back(rname);
								num_valid_read_bases += readp->get_read_length();
								compatible = true;
								break;
							}
						}
					}
					//if (compatible) {
					//	for (unsigned long k = 0; k < read_il.get_num_intervals();
					//			++k) {
					//		ofs << gp->chrom << "\t"
					//			<< "solve\tread_mapping\t"
					//			<< (read_il.get_starts()[k]+1) << "\t"
					//			<< read_il.get_ends()[k] << "\t"
					//			<< ".\t"
					//			<< gp->strand << "\t"
					//			<< ".\t"
					//			<< rname << endl;
					//	}
					//}
				}
				++itr;
			}

			// construct delta_G matrix
			unsigned long num_valid_reads = valid_rnames.size();
			L_(debug) << "Found " << num_valid_reads << " generated read(s)";
			ublas::matrix<double> m_delta_G = ublas::zero_matrix<double>(
					num_valid_reads, isop->num_known_isoforms);
			unsigned long i = 0;
			BOOST_FOREACH(string const & rname, valid_rnames) {
				interval_list<long> const & read_il =
					rname2rinfo_vec[m][rname].get<2>();
				readp->build(read_il, gname2exons[gname].exonps);
				for (unsigned long j = 0;
						j < isop->num_known_isoforms; ++j) {
					if (readp->is_compatible_with_iso(j)) {
						if ((double)readp->get_read_length()
								/ read_il.compute_total_length() > 0.98) {
							m_delta_G(i,j) = readp->prob_generated_by_iso(j);
						}
					}
				}
				++i;
			}
			if (num_valid_reads > 0) {
				m_delta_Gs.push_back(m_delta_G);
			}
			supports.push_back(num_valid_reads);
			support_bases.push_back(num_valid_read_bases);
		}
		//		ofs.close();

		// solve MLE using EM
		vector<double> theta;
		if (m_delta_Gs.size() == 0) {
			theta.resize(isop->num_known_isoforms,
					1.0 / (double)isop->num_known_isoforms);
		} else if (isop->num_known_isoforms == 1) {
			theta.push_back(1);
		} else {
			compute_iso_probs_EM(isop->num_known_isoforms,
					m_delta_Gs, theta);
		}

		// compute RPKM accordingly
		vector<double> rpkm(isop->num_known_isoforms, 0.0);
		double total_read_mbases = 0.0;
		for (unsigned long m = 0; m < M; ++m) {
			total_read_mbases += total_read_bases_vec[m] / 1.0E6;
			for (unsigned int i = 0; i < isop->num_known_isoforms; ++i) {
				rpkm[i] += (double)(support_bases[m]) * theta[i];
			}
		}
		for (unsigned int i = 0; i < isop->num_known_isoforms; ++i) {
			rpkm[i] /= ((double)(
						isop->known_iso_exon_total_lengths[i].back()) / 1.0E3);
			rpkm[i] /= total_read_mbases;
		}

		// compute average log(Pr(s|I, Theta))
		double logll = reads_log_likelihood(isop->num_known_isoforms,
				m_delta_Gs, theta);
		double sum_supports = accumulate(supports.begin(), supports.end(), 0.0);

		L_(debug) << "Estimated Theta and RPKMs for gene " << gname;
		for (unsigned int i = 0; i < isop->num_known_isoforms; ++i) {
			if (theta[i] > 1E-3) {
				L_(debug) << isop->known_iso_names[i] << ": "
					<< theta[i] << ", RPKM: " << rpkm[i]
					<< ", with per read log likelihood: "
					<< ((sum_supports > 1E-5) ? (logll / sum_supports) : 0);
			}
			cout << gname << "\t";
			BOOST_FOREACH(double support, supports) {
				cout << support << "\t";
			}
			cout << isop->known_iso_names[i] << "\t" << theta[i]
				<< "\t" << rpkm[i];
			if (sum_supports > 1E-5) {
				cout << "\t" << (logll / sum_supports) << endl;
			} else {
				cout << "\t" << 0 << endl;
			}
		}

		if (++num_processed_genes > 0 && num_processed_genes % 100 == 0) {
			L_(info) << "Processed " << num_processed_genes << " genes...";
		}
	}
	L_(info) << "Processed " << num_processed_genes << " genes... Done";

	return 0;
}

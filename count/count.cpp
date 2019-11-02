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
	L_(error) << "Usage:\ncount\n	log_level(0,1,2,...) proj_name out_prefix\n	isoform_format isoforms_path g2i_format g2i_path gene_begin_idx gene_end_idx\n  (read_format read_type expected_read_length reads_path )+";
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
		static const string MRF_SINGLE;
};
const string FileFormats::UCSC_GENE2ISOFORM = "UCSC_GENE2ISOFORM";
const string FileFormats::LH_GENE_TXT = "LH_GENE_TXT";
const string FileFormats::MRF_SINGLE = "MRF_SINGLE";

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
	if (argc < 14) {
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
	vector<unsigned long> expected_read_lengths;
	try {
		gene_begin_idx = lexical_cast<unsigned long>(argv[argi++]);
		gene_end_idx = lexical_cast<unsigned long>(argv[argi++]);

		// load read parameters
		while (argi < argc) {
			if (argc - argi < 4) {
				return syntax_err();
			}
			read_formats.push_back(argv[argi++]);
			read_types.push_back(argv[argi++]);
			expected_read_lengths.push_back(
					lexical_cast<unsigned long>(argv[argi++]));
			reads_paths.push_back(argv[argi++]);
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
		if (read_format == FileFormats::MRF_SINGLE) {
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
		} else {
			return unknown_format_err(read_format);
		}
		ifs.close();
		L_(info) << "Sampling method #" << m
		 << ": loaded "	<< rnames_vec[m].size()
		 << " reads associated with the selected gene regions";
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
		vector<double> iso_count(isop->num_known_isoforms, 0.0); //JL
		// record the reads
		//		ofstream ofs;
		//		ofs.clear();
		//		string gene_reads_path = out_prefix + gname + "-reads.gff";
		//		ofs.open(gene_reads_path.c_str());
		//		assert(ofs.is_open() && !ofs.eof());
		//		ofs << "browser position " << isop->chrom
		//			<< ":" << gene_start << "-" << gene_end << endl;
		for (unsigned long m = 0; m < M; ++m) {
		  //			ofs << "track name=\"" << proj_name
		  //				<< ": reads in " << gname
		  //				<< " from sampling method #" << m
		  //				<< "\" visibility=4" << endl;
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
								compatible = true;
								break;
							}
						}
					}
					//					if (compatible) {
					  //						for (unsigned long k = 0; k < read_il.get_num_intervals();
						     //								++k) {
						  //							ofs << gp->chrom << "\t"
							  //								<< "solve\tread_mapping\t"
							  //								<< (read_il.get_starts()[k]+1) << "\t"
							  //								<< read_il.get_ends()[k] << "\t"
							  //								<< ".\t"
							  //								<< gp->strand << "\t"
							  //								<< ".\t"
							  //								<< rname << endl;
							//						}
						//					}
				}
				++itr;
			}

			// count				
			unsigned long num_valid_reads = valid_rnames.size();
			L_(debug) << "Found " << num_valid_reads << " generated read(s)";
			unsigned long i = 0;
			BOOST_FOREACH(string const & rname, valid_rnames) {
				interval_list<long> const & read_il =
					rname2rinfo_vec[m][rname].get<2>();
				readp->build(read_il, gname2exons[gname].exonps);
				for (unsigned long j = 0;
						j < isop->num_known_isoforms; ++j) {
					if (readp->is_compatible_with_iso(j)) {
							iso_count[j]++;     //JL
					}
				}
				++i;
			}
			supports.push_back(num_valid_reads);
		}
		//		ofs.close();

		for (unsigned int i = 0; i < isop->num_known_isoforms; ++i) {
			cout << gname << "\t";
			BOOST_FOREACH(double support, supports) {
				cout << support << "\t";
			}
			cout << isop->known_iso_names[i] << "\t" << iso_count[i] << endl;      //JL
		}

		if (++num_processed_genes > 0 && num_processed_genes % 100 == 0) {
			L_(info) << "Processed " << num_processed_genes << " genes...";
		}
	}
	L_(info) << "Processed " << num_processed_genes << " genes... Done";

	return 0;
}

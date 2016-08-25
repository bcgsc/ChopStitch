/**
 * Find exon exon junctions
 * 
 */

#include "config.h"



#include <utility>
#include <algorithm>
#include <vector>

//#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/adjacency_list.hpp>


//#include "Konnector/konnector.h"
//#include "Bloom/CascadingBloomFilter.h"
//#include "Konnector/DBGBloom.h"
//#include "Konnector/DBGBloomAlgorithms.h"
#include "pstream.h"

//#include "Align/alignGlobal.h"
#include "Common/IOUtil.h"
#include "Common/Options.h"
#include "Common/StringUtil.h"
#include "DataLayer/FastaConcat.h"
#include "DataLayer/FastaInterleave.h"
#include "DataLayer/Options.h"
//#include "Graph/DotIO.h"
//#include "Graph/Options.h"
//#include "Graph/GraphUtil.h"
#include "inotify-cxx.h"

#include "BloomFilter.h"


#include <cassert>
#include <getopt.h>
#include <iostream>
#include <cstring>
#include <algorithm> 
#include <sstream>
#include <fstream>

#if _OPENMP
# include <omp.h>
//# include "Bloom/ConcurrentBloomFilter.h"
#endif

using namespace std;

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define PROGRAM "Expander"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Hamza Khan, Ben Vandervalk\n"

"\n"
"Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k <kmer_size> -o <output_prefix> [options]... <reads1> [reads2]...\n"
"Connect the pairs READS1 and READS2 and close the gap using\n"
"a Bloom filter de Bruijn graph.\n"
"\n"
" Options:\n"
"\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"  -k, --kmer=N               the size of a k-mer\n"
"  -b, --bloom-size=N         size of bloom filter [500M]\n"
"  -X, --internalexons        output internal confident exons in a fasta file \n"
"  -c, --min-coverage=N       kmer coverage threshold for error correction [2].\n"
"                             This option specifies the number of levels in the\n"
"                             cascading Bloom filter; it has no effect if the Bloom\n"
"                             filter is loaded from an external file.\n"
"  -B, --max-branches=N       max branches in de Bruijn graph traversal;\n"
"                             use 'nolimit' for no limit [350]\n"
"  -d, --dot-file=FILE        write graph traversals to a DOT file\n"
"  -D, --dup-bloom-size=N     use an additional Bloom filter to avoid\n"
"                             assembling the same region of the genome\n"
"                             multiple times. This option is highly\n"
"                             recommended when the -E (--extend) option\n"
"                             and has no effect otherwise. As a rule of\n"
"                             thumb, the Bloom filter size should be\n"
"                             about twice the target genome size [disabled]\n"
"  -e, --fixerrors           find and fix single-base errors when reads\n"
"                             have no kmers in bloom filter [disabled]\n"
"  -E, --extend               in addition to finding a connecting path,\n"
"                             extend the reads outwards to the next\n"
"                             dead end or branching point in the de Brujin\n"
"                             graph. If the reads were not successfully\n"
"                             connected, extend them inwards as well.\n"
"  -f, --min-frag=N           min fragment size in base pairs [0]\n"
"  -F, --max-frag=N           max fragment size in base pairs [1000]\n"
"  -i, --input-bloom=FILE     load bloom filter from FILE\n"
"  -I, --interleaved          input reads files are interleaved\n"
"      --mask                 mask new and changed bases as lower case\n"
"      --no-mask              do not mask bases [default]\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -m, --read-mismatches=N    max mismatches between paths and reads; use\n"
"                             'nolimit' for no limit [nolimit]\n"
"  -M, --max-mismatches=N     max mismatches between all alternate paths;\n"
"                             use 'nolimit' for no limit [2]\n"
"  -n  --no-limits            disable all limits; equivalent to\n"
"                             '-B nolimit -m nolimit -M nolimit -P nolimit'\n"
"  -o, --output-prefix=FILE   prefix of output FASTA files [required]\n"
"  -P, --max-paths=N          merge at most N alternate paths; use 'nolimit'\n"
"                             for no limit [2]\n"
"  -q, --trim-quality=N       trim bases from the ends of reads whose\n"
"                             quality is less than the threshold\n"
"      --standard-quality     zero quality is `!' (33)\n"
"                             default for FASTQ and SAM files\n"
"      --illumina-quality     zero quality is `@' (64)\n"
"                             default for qseq and export files\n"
"  -r, --read-name=STR        only process reads with names that contain STR\n"
"  -s, --search-mem=N         mem limit for graph searches; multiply by the\n"
"                             number of threads (-j) to get the total mem used\n"
"                             for graph traversal [500M]\n"
"  -t, --trace-file=FILE      write graph search stats to FILE\n"
"  -v, --verbose              display verbose output\n"
"      --help                 display this help and exit\n"
"      --version              output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

const unsigned g_progressStep = 1000;
/*
 * ignore branches less than this length
 *(false positive branches)
 */
const unsigned g_trimLen = 3;

/*
 * Bloom filter use to keep track of portions
 * of genome that have already been assembled.
 * This Bloom filter is only used when the
 * -E (--extend) option is in effect.
 */
BloomFilter g_dupBloom;

namespace opt {
   
	/** Option for internal exons fasta file **/    
   bool internalexons = false;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** The maximum count value of the Bloom filter. */
	unsigned minCoverage = 2;

	/** Input read files are interleaved? */
	bool interleaved = false;

	/** Max active branches during de Bruijn graph traversal */
	unsigned maxBranches = 350;

	/** multi-graph DOT file containing graph traversals */
	static string dotPath;

	/**
	 * Dup Bloom filter size.
	 * The dup filter is used to avoid assembling duplicate
	 * sequences when the -E (--extend) option is in effect.
	 */
	size_t dupBloomSize = 0;

	/**
	 * Find and fix single base errors when a read has no
	 * kmers in the bloom filter.
	 */
	bool fixErrors = false;

	/**
	 * Extend reads outwards until the next dead or branching
	 * point in the de Bruijn graph.  If a read pair is not
	 * successfully connected, extend them inwards as well.
	 */
	bool extend = false;

	/** The size of a k-mer. */
	unsigned k;

	
	/** Number of hashes. */
	unsigned nhash;
	
	

	/** The minimum fragment size */
	unsigned minFrag = 0;

	/** The maximum fragment size */
	unsigned maxFrag = 1000;

	/** Bloom filter input file */
	static string inputBloomPath;

	/** Max paths between read 1 and read 2 */
	unsigned maxPaths = 2;

	/** Prefix for output files */
	static string outputPrefix;

	/** Max mismatches allowed when building consensus seqs */
	unsigned maxMismatches = 2;

	/** Only process reads that contain this substring. */
	static string readName;

	/** Max mem used per thread during graph traversal */
	static size_t searchMem = 500 * 1024 * 1024;

	/** Output file for graph search stats */
	static string tracefilePath;

	/** Mask bases not in reads */
	static int mask = 0;

	/** Max mismatches between consensus and original reads */
	static unsigned maxReadMismatches = NO_LIMIT;

}

namespace joincount {

   /** Counter for paths found**/
   int foundpath=0;

   /** Counter for no paths found (Failure count)**/
   int nopath=0;

   /** Counter for too many paths found (Failure count)**/
   int toomanypaths=0;

   /** Counter for too many branches found (Failure count)**/
   int toomanybranches=0; 

   /** If the path contains a cycle (Failure count)**/
   int pathscontaincycle=0;

   /** Cases where the search exceeded the assigned memory(Failure count)**/
   int exceedmem = 0; 

}

/** Counters */
static struct {
	size_t noStartOrGoalKmer;
	size_t noPath;
	size_t uniquePath;
	size_t multiplePaths;
	size_t tooManyPaths;
	size_t tooManyBranches;
	size_t tooManyMismatches;
	size_t tooManyReadMismatches;
	size_t containsCycle;
	size_t exceededMemLimit;
	size_t traversalMemExceeded;
	size_t readPairsProcessed;
	size_t readPairsMerged;
	size_t skipped;
	/* counts below are used only when -E is enabled */
	size_t mergedAndExtended;
	size_t mergedAndSkipped;
	size_t singleEndCorrected;
} g_count;

static const char shortopts[] = "b:B:c:d:D:eEf:F:i:Ij:k:lm:M:no:P:q:r:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "min-coverage",     required_argument, NULL, 'c' },
	{ "max-branches",     required_argument, NULL, 'B' },
	{ "dot-file",         required_argument, NULL, 'd' },
	{ "dup-bloom-size",   required_argument, NULL, 'D' },
	{ "fixerrors",       no_argument, NULL, 'e' },
	{ "extend",           no_argument, NULL, 'E' },
	{ "min-frag",         required_argument, NULL, 'f' },
	{ "max-frag",         required_argument, NULL, 'F' },
	{ "input-bloom",      required_argument, NULL, 'i' },
	{ "interleaved",      no_argument, NULL, 'I' },
	{ "threads",          required_argument, NULL, 'j' },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "mask",             no_argument, &opt::mask, 1 },
	{ "no-mask",          no_argument, &opt::mask, 0 },
	{ "no-limits",        no_argument, NULL, 'n' },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "output-prefix",    required_argument, NULL, 'o' },
	{ "read-mismatches",  required_argument, NULL, 'm' },
	{ "max-mismatches",   required_argument, NULL, 'M' },
	{ "max-paths",        required_argument, NULL, 'P' },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "read-name",        required_argument, NULL, 'r' },
	{ "search-mem",       required_argument, NULL, 's' },
	{ "trace-file",       required_argument, NULL, 't' },
	{ "verbose",          no_argument, NULL, 'v' },
   { "internalexons",    no_argument, NULL, 'X' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


/**
  * Global variables
  */
   int big1=0, small1=0, big2=0;
   string read1, read2, final_Seq;
   FastqRecord rec1, rec2, endRead, rec1R, rec1L, rec2R, rec2L;

/**
 * Return true if the Bloom filter contains all of the
 * kmers in the given sequence.
 */
static bool bloomContainsSeq(const BloomFilter& bloom, const Sequence& seq)
{
	if (containsAmbiguityCodes(seq)) {
		Sequence seqCopy = seq;
		flattenAmbiguityCodes(seqCopy, false);
		for (KmerIterator it(seqCopy, opt::k); it != KmerIterator::end();
			++it) {
			if (!bloom[*it])
				return false;
		}
		return true;
	}
	for (KmerIterator it(seq, opt::k); it != KmerIterator::end(); ++it) {
		if (!bloom[*it])
			return false;
	}
	return true;
}

/**
 * Load the kmers of a given sequence into a Bloom filter.
 */
static inline void loadSeq(BloomFilter& bloom, unsigned k, const Sequence& seq)
{
	if (containsAmbiguityCodes(seq)) {
		Sequence seqCopy = seq;
		Sequence rc = reverseComplement(seqCopy);
		flattenAmbiguityCodes(seqCopy, false);
		flattenAmbiguityCodes(rc, false);
		Bloom::loadSeq(bloom, k, seqCopy);
		Bloom::loadSeq(bloom, k, rc);
	} else {
		Bloom::loadSeq(bloom, k, seq);
	}
}

/**
 * Extend a read/pseudoread both left and right until
 * we hit the next dead end or branching point in the
 * de Bruijn graph.
 *
 * @param seq sequence to be extended
 * @param k kmer size
 * @param g de Bruijn graph
 * return true if the read was extended in either
 * (or both) directions, false otherwise
 */
template <typename Graph>
static bool extendRead(Sequence& seq, unsigned k, const Graph& g)
{
	ExtendSeqResult result;
	bool extended = false;

	/*
	 * offset start pos to reduce chance of hitting
	 * a dead end on a false positive kmer
	 */
	const unsigned runLengthHint = 3;
	unsigned startPos = getStartKmerPos(seq, k, FORWARD, g,
		runLengthHint);
	if (startPos != NO_MATCH) {
		assert(startPos <= seq.length() - k);
		result = extendSeq(seq, FORWARD, startPos, k, g,
				NO_LIMIT, g_trimLen, opt::mask);
		if (result == ES_EXTENDED_TO_DEAD_END ||
				result == ES_EXTENDED_TO_BRANCHING_POINT ||
				result == ES_EXTENDED_TO_CYCLE) {
			extended = true;
		}
	}

	startPos = getStartKmerPos(seq, k, REVERSE, g, runLengthHint);
	if (startPos != NO_MATCH) {
		assert(startPos <= seq.length() - k);
		result = extendSeq(seq, REVERSE, startPos, k, g,
				NO_LIMIT, g_trimLen, opt::mask);
		if (result == ES_EXTENDED_TO_DEAD_END ||
				result == ES_EXTENDED_TO_BRANCHING_POINT ||
				result == ES_EXTENDED_TO_CYCLE) {
			extended = true;
		}
	}

	return extended;
}

enum ExtendResult { ER_NOT_EXTENDED, ER_REDUNDANT, ER_EXTENDED };

/**
 * Attempt to extend a merged read (a.k.a. pseudoread)
 * outward to the next branching point or dead end in
 * the de Bruijn graph.
 *
 * @param seq pseudoread to be extended
 * @param k kmer size
 * @param g de Bruijn graph in which to perform extension
 * @return ExtendResult (ER_NOT_EXTENDED, ER_EXTENDED,
 * ER_REDUNDANT)
 */
template <typename Graph>
static inline ExtendResult
extendReadIfNonRedundant(Sequence& seq, unsigned k, const Graph& g)
{
	bool redundant = false;
	if (opt::dupBloomSize > 0) {
		/*
		 * Check to see if the current pseudoread
		 * is contained in a region of the genome
		 * that has already been assembled.
		 */
#pragma omp critical(dupBloom)
		redundant = bloomContainsSeq(g_dupBloom, seq);
		if (redundant)
			return ER_REDUNDANT;
	}
	Sequence origSeq = seq;
	bool extended = extendRead(seq, k, g);
	if (opt::dupBloomSize > 0) {
		/*
		 * mark the extended read as an assembled
		 * region of the genome.
		 */
#pragma omp critical(dupBloom)
		{
			/* must check again to avoid race conditions */
			if (!bloomContainsSeq(g_dupBloom, origSeq))
				loadSeq(g_dupBloom, opt::k, seq);
			else
				redundant = true;
		}
		if (redundant)
			return ER_REDUNDANT;
	}
	assert(!redundant);
	if (extended)
		return ER_EXTENDED;
	else
		return ER_NOT_EXTENDED;
}

/**
 * Print progress stats about reads merged/extended so far.
 */
static inline void printProgressMessage()
{
	cerr << "Merged " << g_count.uniquePath + g_count.multiplePaths << " of "
		<< g_count.readPairsProcessed << " read pairs";

	if (opt::extend) {
		cerr << ", corrected/extended " << g_count.singleEndCorrected << " of "
			<< (g_count.readPairsProcessed - g_count.uniquePath -
				g_count.multiplePaths) * 2
		<< " unmerged reads";
	}

	cerr << " (no start/goal kmer: " << g_count.noStartOrGoalKmer << ", "
		<< "no path: " << g_count.noPath << ", "
		<< "too many paths: " << g_count.tooManyPaths << ", "
		<< "too many branches: " << g_count.tooManyBranches << ", "
		<< "too many path/path mismatches: " << g_count.tooManyMismatches << ", "
		<< "too many path/read mismatches: " << g_count.tooManyReadMismatches << ", "
		<< "contains cycle: " << g_count.containsCycle << ", "
		<< "skipped: " << g_count.skipped
		<< ")\n";
}


/**
 * For a successfully merged read pair, get the sequence
 * representing the connecting path between the two reads.
 */
template <typename Bloom>
static inline string getConnectingSeq(ConnectPairsResult& result,
	unsigned k, const Bloom& bloom)
{
	assert(result.pathResult == FOUND_PATH);
	(void)bloom;

	vector<FastaRecord>& paths = result.mergedSeqs;
	assert(paths.size() > 0);

	Sequence& seq = (paths.size() == 1) ?
		paths.front().seq : result.consensusSeq.seq;

	/*
	 * initialize sequence to the chars between the
	 * start and goal kmers of the path search.
	 */
	int startPos = result.startKmerPos;
	int endPos = seq.length() - result.goalKmerPos - k;
	assert(startPos >= 0 && startPos <=
		(int)(seq.length() - k + 1));

	return seq.substr(startPos, endPos - startPos + k);
}

/** Connect a read pair. */
template <typename Graph, typename Bloom>
static void connectPair(const Graph& g,
	const Bloom& bloom,
	FastqRecord& read1,
	FastqRecord& read2,
	const ConnectPairsParams& params,
	ofstream& mergedStream,
	ofstream& read1Stream,
	ofstream& read2Stream,
	ofstream& traceStream)
{
	/*
	 * Implements the -r option, which is used to only
	 * process a subset of the input read pairs.
	 */
	if (!opt::readName.empty() &&
		read1.id.find(opt::readName) == string::npos) {
#pragma omp atomic
		++g_count.skipped;
		return;
	}

	ConnectPairsResult result =
		connectPairs(opt::k, read1, read2, g, params);

	vector<FastaRecord>& paths = result.mergedSeqs;
	bool mergedSeqRedundant = false;
	bool read1Corrected = false;
	bool read1Redundant = false;
	bool read2Corrected = false;
	bool read2Redundant = false;

	/*
	 * extend reads inwards or outwards up to the
	 * next dead end or branching point in the de
	 * Brujin graph
	 */
	if (opt::extend) {
		ExtendResult extendResult;
		if (result.pathResult == FOUND_PATH
			&& result.pathMismatches <= params.maxPathMismatches
			&& result.readMismatches <= params.maxReadMismatches) {
			assert(paths.size() > 0);
			Sequence& seq = (paths.size() == 1) ?
				paths.front().seq : result.consensusSeq.seq;
			seq = getConnectingSeq(result, opt::k, bloom);
			extendResult = extendReadIfNonRedundant(
				seq, opt::k, g);
			if (extendResult == ER_REDUNDANT) {
#pragma omp atomic
				g_count.mergedAndSkipped++;
				mergedSeqRedundant = true;
			} else if (extendResult == ER_EXTENDED) {
#pragma omp atomic
				g_count.mergedAndExtended++;
			}
		} else {

			/*
			 * read pair could not be merged, so try
			 * to extend each read individually (in
			 * both directions).
			 */

//std::cerr << "correcting " << read1.id << " (read 1)" << std::endl;
			read1Corrected = correctAndExtendSeq(read1.seq,
				opt::k, g, read1.seq.length(), g_trimLen,
				opt::mask);

			if (read1Corrected) {
#pragma omp atomic
				g_count.singleEndCorrected++;
				extendResult = extendReadIfNonRedundant(read1.seq,
					opt::k, g);
				if (extendResult == ER_REDUNDANT)
					read1Redundant = true;
			}

//std::cerr << "correcting " << read2.id << " (read 2)" << std::endl;
			read2Corrected = correctAndExtendSeq(read2.seq,
				opt::k, g, read2.seq.length(), g_trimLen,
				opt::mask);

			if (read2Corrected) {
#pragma omp atomic
				g_count.singleEndCorrected++;
				extendResult = extendReadIfNonRedundant(read2.seq,
					opt::k, g);
				if (extendResult == ER_REDUNDANT)
					read2Redundant = true;
			}
		}
	}

	if (!opt::tracefilePath.empty())
#pragma omp critical(tracefile)
	{
		traceStream << result;
		assert_good(traceStream, opt::tracefilePath);
	}

	switch (result.pathResult) {

		case NO_PATH:
			assert(paths.empty());
			if (result.foundStartKmer && result.foundGoalKmer)
#pragma omp atomic
				++g_count.noPath;
			else {
#pragma omp atomic
				++g_count.noStartOrGoalKmer;
			}
			break;

		case FOUND_PATH:
			assert(!paths.empty());
			if (result.pathMismatches > params.maxPathMismatches ||
					result.readMismatches > params.maxReadMismatches) {
				if (result.pathMismatches > params.maxPathMismatches)
#pragma omp atomic
					++g_count.tooManyMismatches;
				else
					++g_count.tooManyReadMismatches;
				if (opt::extend) {
					if (read1Corrected || read2Corrected)
#pragma omp critical(mergedStream)
					{
						if (read1Corrected && !read1Redundant)
							mergedStream << (FastaRecord)read1;
						if (read2Corrected && !read2Redundant)
							mergedStream << (FastaRecord)read2;
					}
					if (!read1Corrected || !read2Corrected)
#pragma omp critical(readStream)
					{
						if (!read1Corrected)
							read1Stream << (FastaRecord)read1;
						if (!read2Corrected)
							read1Stream << (FastaRecord)read2;
					}
				} else
#pragma omp critical(readStream)
				{
					read1Stream << read1;
					read2Stream << read2;
				}
			}
			else if (paths.size() > 1) {
#pragma omp atomic
				++g_count.multiplePaths;
				if (!mergedSeqRedundant)
#pragma omp critical(mergedStream)
					mergedStream << result.consensusSeq;
			}
			else {
#pragma omp atomic
				++g_count.uniquePath;
				if (!mergedSeqRedundant)
#pragma omp critical(mergedStream)
					mergedStream << paths.front();
			}
			break;

		case TOO_MANY_PATHS:
#pragma omp atomic
			++g_count.tooManyPaths;
			break;

		case TOO_MANY_BRANCHES:
#pragma omp atomic
			++g_count.tooManyBranches;
			break;

		case PATH_CONTAINS_CYCLE:
#pragma omp atomic
			++g_count.containsCycle;
			break;

		case EXCEEDED_MEM_LIMIT:
#pragma omp atomic
			++g_count.exceededMemLimit;
			break;
	}

	if (result.pathResult != FOUND_PATH) {
		if (opt::extend) {
			if (read1Corrected || read2Corrected)
#pragma omp critical(mergedStream)
			{
				if (read1Corrected && !read1Redundant)
					mergedStream << (FastaRecord)read1;
				if (read2Corrected && !read2Redundant)
					mergedStream << (FastaRecord)read2;
			}
			if (!read1Corrected || !read2Corrected)
#pragma omp critical(readStream)
			{
				if (!read1Corrected)
					read1Stream << (FastaRecord)read1;
				if (!read2Corrected)
					read1Stream << (FastaRecord)read2;
			}
		} else
#pragma omp critical(readStream)
		{
			read1Stream << read1;
			read2Stream << read2;
		}
	}
}

/** Connect read pairs. */
template <typename Graph, typename FastaStream, typename Bloom>
static void connectPairs(const Graph& g,
	const Bloom& bloom,
	FastaStream& in,
	const ConnectPairsParams& params,
	ofstream& mergedStream,
	ofstream& read1Stream,
	ofstream& read2Stream,
	ofstream& traceStream)
{
#pragma omp parallel
	for (FastqRecord a, b;;) {
		bool good;
#pragma omp critical(in)
		good = in >> a >> b;
		if (good) {
			connectPair(g, bloom, a, b, params, mergedStream, read1Stream,
				read2Stream, traceStream);
#pragma omp atomic
			g_count.readPairsProcessed++;
			if (opt::verbose >= 2)
#pragma omp critical(cerr)
			{
				if(g_count.readPairsProcessed % g_progressStep == 0)
					printProgressMessage();
			}
		} else {
			break;
		}
	}
}

/**
 * Set the value for a commandline option, using "nolimit"
 * to represent NO_LIMIT.
 */
static inline void setMaxOption(unsigned& arg, istream& in)
{
	string str;
	getline(in, str);
	if (!in.fail() && str.compare("nolimit")==0) {
		arg = NO_LIMIT;
	} else {
		istringstream ss(str);
		ss >> arg;
		// copy state bits (fail, bad, eof) to
		// original stream
		in.clear(ss.rdstate());
	}
}


template <typename Graph, typename Bloom>
static string connectReadPair(const Graph& g,
	const Bloom& bloom,
   FastqRecord& rec1L,
	FastqRecord& read1,
	FastqRecord& read2,
   FastqRecord& rec2L,
   const ConnectPairsParams& params)
{

	ConnectPairsResult result =
		connectPairs(opt::k, read1, read2, g, params);
  
  //cout << result;  

  if(bloom){}

   vector<FastaRecord>& paths = result.mergedSeqs;
   string r1L = rec1L.seq;
   string r1 = read1.seq;
   string r2 = reverseComplement(read2.seq);
   string r2L = reverseComplement(rec2L.seq);
   string return_string;
	switch (result.pathResult) {

		case NO_PATH:
			//assert(paths.empty());
         //cout << "No paths found; insert Ns for "<<read1.seq << "\nNNNNNNNNNN\n"<< read2.seq << endl;   
         //cout << "\nNo paths found\n";
         return_string = r1L + r1 + "NNNNN" + r2 + r2L;
         joincount::nopath+=1;
         //cout << return_string << endl;
			break;

		case FOUND_PATH:
         joincount::foundpath+=1;
			//assert(!paths.empty());
         //cout << "Yuhooo!! These two reads merged \n" << read1 << "\nNNNNNNNNNN\n" << read2 << "\n This is the merged sequence = " << result.consensusSeq << endl;
         //cout << "\nPATH found\n";
			if (paths.size() > 1)
         	return_string = r1L+(string)(result.consensusSeq)+r2L;
			else
				return_string = r1L+(string)(result.mergedSeqs.front().seq)+r2L;
         //cout << return_string<< endl;
			break;

		case TOO_MANY_PATHS:
         //cout << "Too many paths found; inserted Ns for\n "<< read1 << "\nNNNNNNNNNN\n" << read2 << endl;
         //cout << "\nTOO MANY paths found\n";
         joincount::toomanypaths+=1;
         return_string = r1L + r1 + "NNNNN" + r2 + r2L;
         //cout << return_string<< endl;
			break;

		case TOO_MANY_BRANCHES:
         //cout << "Too many branches found; inserted Ns for \n "<< read1 << "\nNNNNNNNNNN\n" << read2 << endl;
         return_string = r1L + r1 + "NNNNN" + r2 + r2L;
        // cout << "\nTOO many Branches\n";
        // cout << return_string<< endl;
         joincount::toomanybranches+=1;
			break;

		case PATH_CONTAINS_CYCLE:
         //cout << "Paths contain cycles; inserted Ns for \n"<< read1 << "\nNNNNNNNNNN\n" << read2 << endl;
         return_string = r1L + r1 + "NNNNN" + r2 + r2L;
        // cout << "\nPaths contain cycles\n";
        // cout << return_string << endl;
         joincount::pathscontaincycle+=1;
			break;

		case EXCEEDED_MEM_LIMIT:
         //cout << "Exceeded memory limit for these two reads; inserted Ns for \n"<< read1 << "\nNNNNNNNNNN\n" <<read2 << endl;
         return_string = r1L + r1 + "NNNNN" + r2 + r2L;
        // cout << "\nExceeded mem_limit\n";
        // cout << return_string <<endl;
         joincount::exceedmem+=1;
			break;
	}
return return_string;

}


void findJunctions(BloomFilter& bloom, int optind, char** argv)

{
   FastaReader reader(argv[optind], FastaReader::FOLD_CASE);
   BloomFilter& bf = bloom;
   
	
	//DBGBloom<BloomFilter> g(*bloom);

	std::ofstream bfile("boundaries.csv", std::ios::app);
   std::ofstream efile("exons.fa", std::ios::app);
   std::ofstream cefile("confident_exons.fa", std::ios::app);


	for (FastaRecord rec; reader >> rec;) 
    {
     if(int((rec.seq).length())<(int(opt::k)+5)){continue;}
	  unsigned pos = (opt::k)-1,
              myflag=0, 
              start = 0, 
              end = 0, 
              n = 10, 
              min_exon=5, 
              unmatch_count=0, 
              snp_chance_A=0,
              snp_chance_T=0,
              snp_chance_G=0,
              snp_chance_C=0,
              snp_nochance=0;

     double FDR = 0.02, 
            leniency = ceil(n*FDR*(opt::k));

	  int first_match=0, 
         last_match_pos=0;

     std::string current_sec = rec.seq;

     //cout << ">" << rec.id << std::endl;
     vector<int> arr;
     vector<int> carr;     

     int str_length = current_sec.length();

     for (unsigned int x=0; x < (str_length-((opt::k)-1)-2); x++) 
            {
             //Kmer it = Kmer(current_sec.substr (x,(opt::k)));
             
		       pos=pos+1;
		       if((bf[it]))
		           {
		               last_match_pos=pos;
		           }
		        
              	
		        if((!bf[it]) && (myflag==0) && (first_match==1))
		           {
		               myflag=1;
		               start = pos;
                     unmatch_count=1;
                     snp_chance_A=0, snp_chance_C=0, snp_chance_G=0, snp_chance_T=0, snp_nochance=0;
                     
		           }
		         if((!bf[it]) && (myflag==1))
		            {
                    if(unmatch_count<=leniency)
                        {  
                         
                         string it_a = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "A");
 
                         string it_t = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "T");
                        
                         string it_g = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "G");
                        
                         string it_c = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "C");
                        

                         if(bf[Kmer(it_a)]){snp_chance_A++;unmatch_count+=1;
                                                 goto SNP_found;}

                         else if(bf[Kmer(it_t)]){snp_chance_T++;unmatch_count+=1;
                                                 goto SNP_found;}
                         
                         else if(bf[Kmer(it_g)]){snp_chance_G++;unmatch_count+=1;
                                                 goto SNP_found;}
                         
                         else if(bf[Kmer(it_c)]){snp_chance_C++;unmatch_count+=1;
                                                 goto SNP_found;}
                         
                         else{snp_nochance++;}
                         unmatch_count+=1;
                       }
                      
                      if(unmatch_count == leniency+1)
                        {

                         if((snp_chance_A  >= 2*snp_nochance) && (snp_chance_A > (snp_chance_T+snp_chance_G+snp_chance_C)))
                             {
                              current_sec.replace(pos-(1+leniency),1,"A");

                             }   
                         if((snp_chance_T  >= 2*snp_nochance) && (snp_chance_T > (snp_chance_A+snp_chance_G+snp_chance_C)))
                             {
                              current_sec.replace(pos-(1+leniency),1,"T");

                             }   
                         if((snp_chance_G  >= 2*snp_nochance) && (snp_chance_G > (snp_chance_A+snp_chance_T+snp_chance_C)))
                             {
                              current_sec.replace(pos-(1+leniency),1,"G");

                              }   
                         if((snp_chance_C  >= 2*snp_nochance) && (snp_chance_C > (snp_chance_A+snp_chance_G+snp_chance_T)))
                             {
                              current_sec.replace(pos-(1+leniency),1,"C");
   
                             }   
                          unmatch_count+=1;
                       }
                                                 
		            }                
                SNP_found:

                if((bf[it]) && (myflag==1) && (first_match==1))
                    {
                       unsigned int next_kmer = x;
                       if (!(bf[Kmer(current_sec.substr (++next_kmer,(opt::k)))]))
	                        {  
                             unmatch_count+=1;
                           }
                    }


              unsigned int a = x;                           

		        if((bf[it]) && (myflag==1) &&  (bf[Kmer(current_sec.substr (++a,(opt::k)))]) && (first_match==1))
	               {  
                       
                        if((bf[Kmer(current_sec.substr (++a,(opt::k)))]))
                            {
                               
			                      myflag = 0;
			                      end = pos;
			                    
                               if(end-start > ((opt::k)+min_exon))
                                    {     
                                          int minus_counter = 0, snp_chance_A=0, snp_chance_C=0, snp_chance_G=0, snp_chance_T=0, snp_nochance=0;
                                          for(int a_minus=(x-1); a_minus>=(x-leniency); a_minus-- )
                                          {
                                           
                                          string it_rev = Kmer(current_sec.substr (a_minus,(opt::k))).str();                               
                                        
                                          string it_a = (it_rev).replace((minus_counter), 1, "A");
                                         
                                          string it_t = (it_rev).replace((minus_counter), 1, "T");
                                      
                                          string it_g = (it_rev).replace((minus_counter), 1, "G");
                                          
                                          string it_c = (it_rev).replace((minus_counter), 1, "C");
                                         
                                          minus_counter = minus_counter + 1;

                                         if(bf[Kmer(it_a)]){snp_chance_A++;}

                                         else if(bf[Kmer(it_t)]){snp_chance_T++;}
                         
                                         else if(bf[Kmer(it_g)]){snp_chance_G++;}
                         
                                         else if(bf[Kmer(it_c)]){snp_chance_C++;}
                         
                                         else{snp_nochance++;}

                                         }

                                         if(snp_chance_A  >= 2*snp_nochance || snp_chance_T  >= 2*snp_nochance || snp_chance_G  >= 2*snp_nochance || snp_chance_C  >= 2*snp_nochance )
                                             {
                                          
                                           bfile << "," << start-1 << "\n" << (rec.id) << ","<< start ;
                                           arr.push_back(start-1);
                                           arr.push_back(start);    
                                             }
                                         else 
                                             {
                                       
                                           bfile << "," << start-1 << "\n" << (rec.id) <<","<< start << "," << (end-opt::k)-1<< "\n" << (rec.id) <<","<< end-opt::k ;
                                           arr.push_back(start-1);
                                           arr.push_back(start);
                                           arr.push_back((end-opt::k)-1);
                                           arr.push_back(end-opt::k);
                                             }     
                                        }
                    
			                      if((end - start >= (((opt::k)-leniency)-2)) && ((end-start)<= ((opt::k)+min_exon)))
                                     {
			                                 
			                                  bfile << "," << start-1 << "\n" << (rec.id) <<","<< end-opt::k ;
                                           arr.push_back(start-1);
                                           arr.push_back(end-opt::k);
			                            }
                             }
		             }

              //Testing the first matching kmer in the bloom filter
              if((bf[it]) && first_match==0)
		              {
                          unsigned int next_it = x;
                          //Added extra checks to avoid False positives
                          if(bf[Kmer(current_sec.substr (++next_it,(opt::k)))] && bf[Kmer(current_sec.substr (++next_it,(opt::k)))])
                            { 
                             bfile << (rec.id)<< ",";
	    		                 bfile << pos-((opt::k)-1);
                             arr.push_back(pos-((opt::k)-1));
			                    first_match=1;
                            }
		              }
        
		    }  
      if(first_match==1)
         {
		      	bfile << "," << last_match_pos << "\n";
               arr.push_back(last_match_pos);
         } 

     

  if(opt::internalexons)
         {
             int check_flag=0;
             carr=arr;
            /*
             for (std::vector<int>::const_iterator i = carr.begin(); i != carr.end(); ++i)
             std::cout << *i << ',';
             for (std::vector<int>::const_iterator i = arr.begin(); i != arr.end(); ++i)
             std::cout << *i << ',';
            */

             if(!carr.empty() && carr.size()==2)
                 {
                   if (*(carr.begin())!=1 && *(carr.end()-1) < str_length-2)
                      {
                      string single_read = rec.seq.substr(((*(carr.begin()))-1),((*(carr.end()-1)-(*(carr.begin()))+1)));
                      cefile << ">" << rec.id << "_" << *(carr.begin())<<"_"<< *(carr.end()-1)<< std::endl;
                      cefile << single_read<<endl; 
                      std::flush(cefile);             
                      }              
                  }

            if(!carr.empty() && (carr.size())>2) 
                 {   
                  int beGIN = *(carr.begin());
                  int enD = *(carr.end()-1);
                  check_flag=1;

                  if(beGIN == 1)
                     {
                       carr.erase(carr.begin(),carr.begin()+2);
                     }
                 
                  if(enD >= str_length-2)
                     {
                       
                       carr.erase (carr.end()-2,carr.end());
                     }
                  }

           
             if(!carr.empty() && carr.size()==2 && check_flag==1)
                 {

                      string single_read = rec.seq.substr(((*(carr.begin()))-1),((*(carr.end()-1)-(*(carr.begin()))+1)));
                      cefile << ">" << rec.id << "_" << *(carr.begin())<<"_"<< *(carr.end()-1)<< std::endl;
                      cefile << single_read<<endl; 
                      std::flush(cefile);             
                      check_flag=0;             
                  }


            if(!carr.empty() && (carr.size())>2) 
                 { 
                     int beGIN = *(carr.begin());
                     int enD = *(carr.end()-1);
                     int Big1=0, Small1=0, Big2=0;
                     string Read1, Read2;

                     for (std::vector<int>::const_iterator w = (carr.begin())+1; w != (carr.end())-1; w+=2)
                         {
                         if(*w < *(w+1))
                             {
                                  Big1 = *(w+1);
                                  Small1 = *w;
                             }
                         else
                             {
                                  Small1 = *(w+1);
                                  Big1 = *w;   
                             }
              
                         if(w == (carr.begin())+1) 
                             {
                                 Read1 = rec.seq.substr(beGIN-1,(Small1-beGIN)+1);
                                 cefile << ">" << rec.id << "_" << beGIN<<"_"<< Small1<< std::endl;
                                 cefile << Read1<<endl;
                                 std::flush(cefile); 
                                 Big2 = Big1; 
                             } 
                         else 
                             {
                                 Read2 = rec.seq.substr(Big2-1,(Small1-Big2)+1);
                                 cefile << ">" << rec.id << "_" << Big2<<"_"<< Small1<< std::endl;
                                 cefile << Read2<<endl;
                                 std::flush(cefile); 
                                 Big2 = Big1;  
                             }       
                         }
        
                         if(!carr.empty() && (carr.size())>2) 
                             {                    
                                 string End_Read = rec.seq.substr(Big2-1,(enD-Big2)+1);
                                 cefile << ">" << rec.id << "_" << Big2<<"_"<< enD<< std::endl;
                                 cefile << End_Read<<endl;
                                 std::flush(cefile);                   
                             } 
                  }


         }


    if(!arr.empty() && arr.size()==2)
        {
             
               string single_read = rec.seq.substr(((*(arr.begin()))-1),((*(arr.end()-1)-(*(arr.begin()))+1)));
               efile << ">" << rec.id << "_" << *(arr.begin())<<"_"<< *(arr.end()-1)<< std::endl;
               efile << single_read<<endl; 
               std::flush(efile);             
              
        }


     if(!arr.empty() && (arr.size())>2) 
        {   
              int bEgin = *(arr.begin());
              int eNd = *(arr.end()-1);

              for (std::vector<int>::const_iterator q = (arr.begin())+1; q != (arr.end())-1; q+=2)
                 {
                   if(*q < *(q+1))
                      {
                        big1 = *(q+1);
                        small1 = *q;
 
                      }
                  else
                      {
                         small1 = *(q+1);
                         big1 = *q;
                        
                      }
              
                  if(q == (arr.begin())+1) 
                      {
                           read1 = rec.seq.substr(bEgin-1,(small1-bEgin)+1);
                           efile << ">" << rec.id << "_" << bEgin<<"_"<< small1<< std::endl;
                           efile << read1<<endl;
                           std::flush(efile); 
                           big2 = big1; 
                           
                      } 
                 else 
                      {
                          
                            read2 = rec.seq.substr(big2-1,(small1-big2)+1);
                            efile << ">" << rec.id << "_" << big2<<"_"<< small1<< std::endl;
                            efile << read2<<endl;
                            std::flush(efile); 
                            big2 = big1;  

                   } 

              
           }
        
         if(!arr.empty() && (arr.size())>2) 
               {                    
                     string end_read = rec.seq.substr(big2-1,(eNd-big2)+1);
                     efile << ">" << rec.id << "_" << big2<<"_"<< eNd<< std::endl;
                     efile << end_read<<endl;
                     std::flush(efile); 
                   
               } 


  }
}

 efile.close();     


}	


/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */
int main(int argc, char** argv)
{
	bool die = false;
	bool minCovOptUsed = false;

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true; break;
		  case 'b':
			opt::bloomSize = SIToBytes(arg); break;
		  case 'c':
			arg >> opt::minCoverage;
			minCovOptUsed = true;
			break;
		  case 'B':
			setMaxOption(opt::maxBranches, arg); break;
		  case 'd':
			arg >> opt::dotPath; break;
		  case 'D':
			opt::dupBloomSize = SIToBytes(arg); break;
		  case 'e':
			opt::fixErrors = true; break;
		  case 'E':
			opt::extend = true; break;
		  case 'f':
			arg >> opt::minFrag; break;
		  case 'F':
			arg >> opt::maxFrag; break;
		  case 'i':
			arg >> opt::inputBloomPath; break;
		  case 'I':
			opt::interleaved = true; break;
		  case 'j':
			arg >> opt::threads; break;
		  case 'k':
			arg >> opt::k; break;
		  case 'm':
			setMaxOption(opt::maxReadMismatches, arg); break;
		  case 'n':
			opt::maxBranches = NO_LIMIT;
			opt::maxReadMismatches = NO_LIMIT;
			opt::maxMismatches = NO_LIMIT;
			opt::maxPaths = NO_LIMIT;
			break;
		  case 'M':
			setMaxOption(opt::maxMismatches, arg); break;
		  case 'o':
			arg >> opt::outputPrefix; break;
		  case 'P':
			setMaxOption(opt::maxPaths, arg); break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 'r':
			arg >> opt::readName; break;
		  case 's':
			opt::searchMem = SIToBytes(arg); break;
		  case 't':
			arg >> opt::tracefilePath; break;
		  case 'v':
			opt::verbose++; break;
		  case 'X':
			opt::internalexons = true; break;
		  case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
		die = true;
	}

	if (opt::outputPrefix.empty()) {
		cerr << PROGRAM ": missing mandatory option `-o'\n";
		die = true;
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing input file arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (!opt::inputBloomPath.empty() && minCovOptUsed) {
		cerr << PROGRAM ": warning: -c option has no effect when "
			" using a pre-built Bloom filter (-i option)\n";
	}


#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	Kmer::setLength(opt::k);

#if USESEQAN
	seqanTests();
#endif

	assert(opt::bloomSize > 0);

	if (opt::dupBloomSize > 0)
		g_dupBloom.resize(opt::dupBloomSize * 8);

	//BloomFilter* bloom;
	//CascadingBloomFilter* cascadingBloom = NULL;
	BloomFilter bloom;
	

	if (!opt::inputBloomPath.empty()) {

		if (opt::verbose)
			std::cerr << "Loading Bloom filter from `"
				<< opt::inputBloomPath << "'...\n";

		//bloom = new BloomFilter();

		const char* inputPath = opt::inputBloomPath.c_str();
		//ifstream inputBloom(inputPath, ios_base::in | ios_base::binary);
		//assert_good(inputBloom, inputPath);
		//inputBloom >> *bloom;
		//assert_good(inputBloom, inputPath);
		//inputBloom.close();
		bloom = new BloomFilter(sbfSize*opt::ibits, opt::nhash, opt::k,inputPath);
		

	} else {

		if (opt::verbose)
			std::cerr << "Using a minimum kmer coverage threshold of "
				<< opt::minCoverage << "\n";

		// Specify bloom filter size in bits and divide by number
		// of levels in cascading Bloom filter.

		size_t bits = opt::bloomSize * 8 / opt::minCoverage;
		cascadingBloom = new CascadingBloomFilter(bits, opt::minCoverage);
#ifdef _OPENMP
		ConcurrentBloomFilter<CascadingBloomFilter> cbf(*cascadingBloom, 1000);
		for (int i = optind; i < argc; i++)
			Bloom::loadFile(cbf, opt::k, string(argv[i]), opt::verbose);
#else
		for (int i = optind; i < argc; i++)
			Bloom::loadFile(*cascadingBloom, opt::k, string(argv[i]), opt::verbose);
#endif
		bloom = &cascadingBloom->getBloomFilter(opt::minCoverage - 1);
	}
	
	if (opt::verbose)
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< 100 * bloom->FPR() << "%\n";


	FastaReader freader(argv[optind], FastaReader::FOLD_CASE);

   findJunctions(*bloom,optind, argv);



	if (!opt::inputBloomPath.empty())
		delete bloom;
	else
		delete cascadingBloom;

	return 0;
}


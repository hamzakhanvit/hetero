#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <cassert>
#include <cmath>
#include <utility>
#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <cstring>
#include "BloomFilter.hpp"
#include "pstream.h"
#include "Options.h"
#include "FastaConcat.h"
#include <unordered_map>

#define PROGRAM "hetero"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamza Khan\n"
    "Copyright 2018 OICR\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Cluster haploid reads based on kmer fasta file\n"
    "Accepatble file formats: fasta\n"
    "\n"
    " Options:\n"
    "\n"
    "  -i, --input-bloom=FILE     load bloom filter from FILE\n" 
    "      --help	          display this help and exit\n"
    "      --version	          output version information and exit\n"
    "\n"
    "Report bugs to Hamza.Khan@oicr.on.ca\n";

using namespace std;

namespace opt {
unsigned k=50;
unsigned nhash;
static string inputBloomPath;
size_t m;
double FPR=0.0;
}

static const char shortopts[] = "i:hv";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "input-bloom",	required_argument, NULL, 'i' },  
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


void kmertoreads(BloomFilter& bloom, int optind, char** argv) {
    FastaReader reader(argv[optind], FastaReader::FOLD_CASE);
    BloomFilter& bf = bloom;
    std::ofstream kfile("kmersinreads.tsv");
    kfile << "kmer\tReadIDs";
      
    unordered_map<string, vector<string>> kmerreads;

    for (FastaRecord rec; reader >> rec;) {
        if(int((rec.seq).length())<(int(opt::k)+5))
            continue;

        std::string current_sec = rec.seq;
        unsigned int pos = (opt::k)-1;
        int str_length = current_sec.length();

        for (unsigned int x=0; x < (str_length-((opt::k)-1)-2); x++) {
            string it = current_sec.substr(x,(opt::k));
            pos=pos+1;
            if(bf.contains(it.c_str()))
                std::cout << "Present kmer " << it << " in "  << rec.id << " Barcode = "<< rec.comment <<"\n";
                kmerreads[it].push_back(rec.id+":"+rec.comment);
        }
     }
    
      for (auto it = kmerreads.begin(); it != kmerreads.end(); it++)
         {
         kfile << it->first << "\t";
         for(std::size_t it2 = 0; it2 < (*it).second.size(); ++it2)
             kfile <<(*it).second[it2] << ",";
         kfile <<"\n";
     } 

}


int main(int argc, char** argv) {

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 'i':
            arg >> opt::inputBloomPath;
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (argc - optind < 0) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    if (!opt::inputBloomPath.empty()) {

        if (opt::verbose)
            std::cerr << "Loading Bloom filter from `"
                      << opt::inputBloomPath << "'...\n";

        string infoPath;
        size_t pos = opt::inputBloomPath.rfind("/");
        if (pos == std::string::npos)
            infoPath = "Bfilter.inf";
        else {
            infoPath=opt::inputBloomPath.substr(0, pos) + "/Bfilter.inf";
        }

        cerr << infoPath << "\n";

        ifstream bfinfo(infoPath.c_str());
        string bfLine;
        while(bfinfo >> opt::m >> opt::nhash >> opt::k >> opt::FPR) {}
        cerr << opt::m << "\t" << opt::nhash<< "\t" <<opt::k <<"\t" << opt::FPR << "\n";
        bfinfo.close();

        BloomFilter bloom(opt::m, opt::nhash, opt::k, opt::inputBloomPath.c_str());
        cerr << bloom.getPop() << "\n";

        kmertoreads(bloom, optind, argv);
    }
    else
        cerr << "No input file\n";
}


#ifdef _OPENMP
# include <omp.h>
#endif

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
#include <typeinfo>
#include <unordered_set>

#define PROGRAM "hetero"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamza Khan\n"
    "Copyright 2018 OICR\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Cluster haploid reads based on kmer fasta file\n"
    "Acceptable file formats: fasta, fastq\n"
    "\n"
    " Options:\n"
    "\n"
    "  -i, --input-bloom=FILE     load bloom filter from FILE\n"
    "  -f, --kmerfasta=FILE       kmer fasta FILE\n" 
    "-t, --threads=N                    Use N parallel threads [1]\n" 
    "      --help	          display this help and exit\n"
    "      --version	          output version information and exit\n"
    "\n"
    "Report bugs to Hamza.Khan@oicr.on.ca\n";

using namespace std;

namespace opt {
unsigned k=50;
unsigned nhash;
unsigned nThrd=1;
static string inputBloomPath;
static string kmerfastafile;
size_t m;
double FPR=0.0;
}

namespace hetero {
unordered_set<string> hetero_kmers;

}

static const char shortopts[] = "t:i:f:hv";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "input-bloom", required_argument, NULL, 'i' },  
    { "kmerfasta", required_argument, NULL, 'f' },
    {"threads", required_argument, NULL, 't' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


bool is_kmer_hetero(string it, unordered_set<string>& kmers){
//bool is_kmer_hetero(string it, BloomFilter& bf){
    
    int hetero_chance = 0;
    string temp;
    //cout << "kmer =" << it << endl;

    vector<char> bases = {'A','T','G','C'};
    //vector<string> check; 
    for (size_t x=0; x < it.length(); x++) {
      
        for (std::vector<char>::iterator i=bases.begin(); i!=bases.end(); i++){
           string c(1, *i);

           temp = it;
           string new_kmer = ((temp).replace(x, 1, c)).c_str();
 
           //if((*i)!=(temp[x])) 
           //if((*i)!=(it.at(x)))
           //   cout << "++++++++++Replaced at " << x << " with " << c << " to give " << new_kmer << endl;
   
           if(((*i)!=(it.at(x))) && (kmers.find(new_kmer) != kmers.end())){
           //if( ((*i)!=(temp[x])) && (bf.contains(((temp).replace(x, 1, c)).c_str()))){  
                cout << "HETERO Replaced at " << x << " with " << c << " to give " << new_kmer << endl;
                //check.push_back(new_kmer);
                hetero_chance++;
            } 
       }
    }

    //for(auto &i:check)
    //    cerr << "\nCheck=" << i << endl;

    //cout << "hetero_chance = " << hetero_chance << endl;
    return hetero_chance == 1 ? 1 : 0;
}


void kmertoreads(BloomFilter& bloom, int optind, char** argv) {
    BloomFilter& bf = bloom;
    std::ofstream kfile("kmersinreads.tsv");
    kfile << "kmer\tReadIDs";
    
    cout  << argv[optind] << endl; 
    FastaReader kmerreader((opt::kmerfastafile).c_str(), FastaReader::FOLD_CASE);
    //FastaReader kmerreader(argv[optind], FastaReader::FOLD_CASE);
    cout << "kmer Filename = " << opt::kmerfastafile << endl;   

    unordered_map<string, vector<string>> kmerreads;

    unordered_set<string> kmers; 

    for (FastaRecord kmer; kmerreader >> kmer;) {
        kmers.insert(kmer.seq);
        //bool is_hetero = is_kmer_hetero(kmer.seq,bf); 
        //cout << kmer.seq << " is hetero? " << is_hetero << endl;
    }

    cerr << "Size of kmers set = " << kmers.size() << endl;

    
    #pragma omp parallel
    {   
        #pragma omp single 
        {
            for (unordered_set<string> :: iterator itr = kmers.begin(); itr != kmers.end(); itr++) {
                #pragma omp task 
                {
                    bool is_hetero = is_kmer_hetero((*itr),kmers);
                    cout << (*itr) << " is hetero? " << is_hetero << endl;
                    if(is_hetero) {
                        hetero::hetero_kmers.insert(*itr);
                    }
                }
            }
        }
            
    }
    

    for (auto &y:hetero::hetero_kmers)
         cerr << "-hetero = " << y << endl;
   
/*
    FastaReader reader(argv[optind], FastaReader::FOLD_CASE);
    cout << "Filename = " << argv[optind] << " typeid(variable).name() = " << typeid(argv[optind]).name()<< endl;
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
                //kmerreads[it].push_back(rec.id+":"+rec.comment);
        }
     }
    
      for (auto it = kmerreads.begin(); it != kmerreads.end(); it++)
         {
         kfile << it->first << "\t";
         for(std::size_t it2 = 0; it2 < (*it).second.size(); ++it2)
             kfile <<(*it).second[it2] << ",";
         kfile <<"\n";
     } 

 */
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
        case 'f':
            arg >> opt::kmerfastafile;
            break;
        case 't':
            arg >> opt::nThrd;
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

    #ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
    #endif

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

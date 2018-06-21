
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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

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
unordered_map<string, vector<string>> kmerbarcodes;
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
                //cout << "HETERO Replaced at " << x << " with " << c << " to give " << new_kmer << endl;
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


void findheteroreads(BloomFilter& bloom, int optind, char** argv) {
    BloomFilter& bf = bloom;
    
    cout  << argv[optind] << endl; 
    FastaReader kmerreader((opt::kmerfastafile).c_str(), FastaReader::FOLD_CASE);
    //FastaReader kmerreader(argv[optind], FastaReader::FOLD_CASE);
    cout << "kmer Filename = " << opt::kmerfastafile << endl;   

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
                    //cout << (*itr) << " is hetero? " << is_hetero << endl;
                    if(is_hetero) {
                        hetero::hetero_kmers.insert(*itr);
                    }
                }
            }
        }
            
    }
    //for (auto &y:hetero::hetero_kmers)
         //cerr << "-hetero = " << y << endl;  

    kmers.clear();
}


void findcommonbarcodes(int optind, char** argv) {  

    std::ofstream kfile("kmersandbarcodes.tsv");
    kfile << "kmer\tbarcodes\n";

    FastaReader reader(argv[optind], FastaReader::FOLD_CASE);
   
    for (FastaRecord rec; reader >> rec;) {
        if(int((rec.seq).length())<(int(opt::k)+5))
            continue;

        std::string current_sec = rec.seq;
        unsigned int pos = (opt::k)-1;
        int str_length = current_sec.length();

        #pragma omp parallel for schedule(guided) shared (hetero::hetero_kmers, hetero::kmerbarcodes)
        for (unsigned int x=0; x < (str_length-((opt::k)-1)-2); x++) {
            string it = current_sec.substr(x,(opt::k));
            pos=pos+1;
            if((hetero::hetero_kmers.find(it.c_str()) != hetero::hetero_kmers.end()) && (rec.comment!="")) {
            //if(bf.contains(it.c_str())) {
                //std::cout << "Present kmer " << it << " in "  << rec.id << " Barcode = "<< rec.comment <<"\n";
                hetero::kmerbarcodes[it].push_back(rec.comment);
            }
        }
    }
    
    for (auto it = hetero::kmerbarcodes.begin(); it != hetero::kmerbarcodes.end(); it++) {
        kfile << it->first << "\t";
        for(std::size_t it2 = 0; it2 < (*it).second.size(); ++it2)
            kfile << (*it).second[it2] << "\t";
        kfile <<"\n";
    } 

}




struct Vertex {
    int id;
    const char* name;
    Vertex(int i = -1, const char* name = "default") : id(i), name(name) {}
};

template <typename It> boost::iterator_range<It> mir(std::pair<It, It> const& p) {
    return boost::make_iterator_range(p.first, p.second);
}

template <typename It> boost::iterator_range<It> mir(It b, It e) {
    return boost::make_iterator_range(b, e);
}

typedef typename boost::adjacency_list<
    boost::setS, boost::vecS,
    boost::undirectedS,
    Vertex,                                      // bundled properties (id, name)
    boost::property<boost::edge_weight_t, float> // interior property
    > Graph;

typedef Graph::edge_descriptor Edge;

void save(Graph const& g, const char* fname);  


#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/format.hpp>
#include <fstream>
void save(Graph const& g, const char* fname) {
    std::ofstream ofs(fname);
    using namespace boost;

    write_graphviz(
            ofs,
            g,
            make_label_writer(make_function_property_map<Graph::vertex_descriptor, std::string>([&] (Graph::vertex_descriptor v){ return g[v].name; })),
            make_label_writer(make_transform_value_property_map([](float v){return boost::format("%1.1f") % v;}, boost::get(edge_weight, g)))
        );
}


void numcomponents(Graph &graph){
    
    std::vector<int> component (boost::num_vertices (graph));
    size_t num_components = boost::connected_components (graph, &component[0]);
    std::cout << "num_components = " << num_components <<std::endl;
    
    /*
    std::cout << "Vertices in the each component:" << std::endl;
    for (size_t i = 0; i < num_components; ++i) { 
        for (size_t j = 0; j < boost::num_vertices (graph); ++j) { 
            if (component[j] == i)
                std::cout << j << " ";     
         }
    }
    */
}


void makegraph(){
    Graph graph;
    
    int vertex_id=-1;
    unordered_map<string,int> barcodeID;
    //boost::graph_traits<Graph>::vertex_descriptor vertexA = vertex(0,graph);
    //boost::graph_traits<Graph>::vertex_descriptor vertexB = vertex(0,graph);
    boost::graph_traits<Graph>::vertex_descriptor vertexA = add_vertex(Vertex {0, "DummyA" }, graph);
    boost::graph_traits<Graph>::vertex_descriptor vertexB = add_vertex(Vertex {1, "DummyB"}, graph);
    //remove_vertex(vertexA, graph);
    //remove_vertex(vertexB, graph);
    graph.clear();


    for (auto it = hetero::kmerbarcodes.begin(); it != hetero::kmerbarcodes.end(); it++) {

        set<string> s( (*it).second.begin(), (*it).second.end() );
        (*it).second.assign( s.begin(), s.end() );

        if((*it).second.size()>1) {

        for(std::size_t it2 = 0; it2 < (*it).second.size(); ++it2){ 
             if (barcodeID.find(((*it).second[it2]).c_str()) == barcodeID.end()) {
                 vertex_id++;
                 vertexA = add_vertex(Vertex {vertex_id, (((*it).second[it2])).c_str() }, graph);
                 barcodeID[((*it).second[it2]).c_str()] = vertex_id;
             }

             for(std::size_t it3 = 0; it3 < (*it).second.size(); ++it3){
                  if(it2!=it3){
                     if (barcodeID.find(((*it).second[it3]).c_str()) == barcodeID.end()) {
                         vertex_id++;
                         vertexB = add_vertex(Vertex { vertex_id,(((*it).second[it3])).c_str() }, graph);
                         barcodeID[((*it).second[it3]).c_str()] = vertex_id;
                      
                      }          
                             
                      if((boost::edge(barcodeID[((*it).second[it2]).c_str()], barcodeID[((*it).second[it3]).c_str()], graph).second) 
                         && (((*it).second[it2])!=((*it).second[it3]))){ 
                          std::pair<Edge, bool> ed = boost::edge(barcodeID[((*it).second[it2]).c_str()], barcodeID[((*it).second[it3]).c_str()],graph);
                          int weight = get(boost::edge_weight_t(), graph, ed.first);
                          //cerr << "Weight = " << weight << endl;
                          int weightToAdd = 1.0;
                          boost::put(boost::edge_weight_t(), graph, ed.first, weight+weightToAdd);
                      }
                      if (!(boost::edge(barcodeID[((*it).second[it2]).c_str()], barcodeID[((*it).second[it3]).c_str()], graph).second) 
                           && (((*it).second[it2])!=((*it).second[it3]))) {
                           //cout << "New edge\n";
                           add_edge(vertexA, vertexB,  1.0f, graph);
                      }
                  }
              }
         }

         }
    }




      
    std::cout << "# of vertices : " << num_vertices(graph) << "\n";
    std::cout << "# of edges:    " << num_edges(graph)    << "\n";
    /*   
    for (auto vd : boost::make_iterator_range(vertices(graph))) {
    std::cout << "Vertex descriptor #" << vd 
         << " degree:" << degree(vd, graph)
         << " id:"     << graph[vd].id
         << " name:"     << graph[vd].name
         << "\n";
    }
    */
    save(graph, "before.dot");
    numcomponents(graph);
  
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

        findheteroreads(bloom, optind, argv);
        findcommonbarcodes(optind, argv);
        makegraph();
    }
    else
        cerr << "No input file\n";
}

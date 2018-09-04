
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
#include "Sequence.h"
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
    "  -t, --threads=N            Use N parallel threads [1]\n"
    "  -k, --kmersize=N           kmer size[41]\n"
    "      --help	          display this help and exit\n"
    "      --version	          output version information and exit\n"
    "\n"
    "Report bugs to Hamza.Khan@oicr.on.ca\n";

using namespace std;

namespace opt {
unsigned k=41;
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

static const char shortopts[] = "t:i:f:k:hv";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "input-bloom", required_argument, NULL, 'i' },  
    { "kmerfasta", required_argument, NULL, 'f' },
    { "kmersize", required_argument, NULL, 'k'},
    { "threads", required_argument, NULL, 't' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


void is_kmer_hetero(string it, unordered_set<string>& kmers){
//bool is_kmer_hetero(string it, BloomFilter& bf){
    
    int hetero_chance = 0;
    //cout << "kmer =" << it << endl;
    std::string test_kmer = it;
    //vector<string> check;

    //vector<char> bases = {'A','T','G','C'};
    const char* bases = "ACGT";

    for (size_t x=0; x < it.length(); x++) {
        ////int FID = omp_get_thread_num();
        ////cout << "x=" << x << " FID = " << FID << endl;     
         
        test_kmer = it;
        char curr_base = test_kmer[x];

        for(size_t bi = 0; bi < 4; ++bi) {
            char test_base = bases[bi];
            if(test_base == curr_base) {
                continue;
            }

            test_kmer[x] = test_base;
            //cout << "tmer =" << test_kmer << endl;
            hetero_chance += kmers.find(canonical_representation(test_kmer)) != kmers.end();
       }
    }

    //for(auto &i:check)
    //    cerr << "\nCheck=" << i << endl;

    ////cout << "hetero_chance = " << hetero_chance << endl;
    int is_hetero = (hetero_chance == 1 ? 1 : 0);
    ////cout << it << " is hetero? " << is_hetero << endl;
    if(is_hetero) {
        #pragma omp critical
        hetero::hetero_kmers.insert(it);
    }
}


void findheteroreads(int optind, char** argv) {
//void findheteroreads(BloomFilter& bloom, int optind, char** argv) {
    //BloomFilter& bf = bloom;
    
    cout  << argv[optind] << endl; 
    unordered_set<string> kmers;
   
    FastaReader kmerreader((opt::kmerfastafile).c_str(), FastaReader::FOLD_CASE);
    cout << "kmer Filename = " << opt::kmerfastafile << endl;    
 
    for (FastaRecord kmer;kmerreader >> kmer;){
        kmers.insert(canonical_representation(kmer.seq));
        //bool is_hetero = is_kmer_hetero(kmer.seq,bf); 
        //cout << kmer.seq << " is hetero? " << is_hetero << endl;
    }
    
    cerr << "Size of kmers set = " << kmers.size() << endl;
    //vector<string> kmers_vec;
    //kmers_vec.insert(kmers_vec.end(), kmers.begin(), kmers.end());

   /* 
    #pragma omp parallel
    {   
        #pragma omp single 
        {
            for (unordered_set<string> :: iterator itr = kmers.begin(); itr != kmers.end(); itr++) {
                #pragma omp task shared (hetero::hetero_kmers, kmers)
                {
                    is_kmer_hetero((*itr),kmers);
                    //cout << (*itr) << " is hetero? " << is_hetero << endl;
                    //if(is_hetero) {
                    //    hetero::hetero_kmers.insert(*itr);
                    //}
                }
            }
        }
            
    }

    */

    /*
    int kmersize = kmers.size();
    #pragma omp parallel for schedule(guided) shared (hetero::hetero_kmers, kmers)
    for(int b = 0; b < kmersize; b++){
        ////int ID = omp_get_thread_num();
        ////cout << "b=" << b << " ID = " << ID << endl;
        unordered_set<string> :: iterator itr = kmers.begin();
        advance(itr, b);
        is_kmer_hetero((*itr),kmers);
        ////cout << "LOOP=" << b << endl;
    }
    */

    //Convert set to vector for openmp parallelization
    vector<string> kmers_vec;
    //kmers_vec.assign( kmers.begin(), kmers.end() );   
    kmers_vec.insert(kmers_vec.end(), kmers.begin(), kmers.end());    
   
     
    int kmers_vec_size = kmers_vec.size();   
    cout << "kmers_vec_size = " << kmers_vec_size << endl;
 
    #pragma omp parallel for schedule(guided) shared (hetero::hetero_kmers, kmers) 
    for (int i = 0; i < kmers_vec_size; i++) {
        //cerr << i << endl;
        is_kmer_hetero(kmers_vec[i],kmers);

    }

    cerr << "Completed findheteroreads" << ", Size of hetero::hetero_kmers = "<<  hetero::hetero_kmers.size()<< endl;
    for (auto &y:hetero::hetero_kmers){
        cerr << "-hetero = " << y << endl;  
    }
    //kmers.clear();
}



void findcommonbarcodes(int optind, char** argv) {  

    std::ofstream kfile("kmersandbarcodes.tsv");
    kfile << "kmer\tbarcodes\n";

    FastaReader reader(argv[optind], FastaReader::FOLD_CASE);
  
    std::ofstream headerfile("headerfile.tsv");
    headerfile << "ID\tbarcode\n";
 
    for (FastaRecord rec; reader >> rec;) {
        if((int)rec.seq.length()<(int)(opt::k+5))
            continue;
        std::string current_sec = rec.seq;
        headerfile<<rec.id<<"\t"<<rec.comment<<endl;
        //cout << rec.id << " - "  << rec.comment <<endl;
        //cout << rec.seq << endl; 
        //unsigned int pos = (opt::k)-1;
        int str_length = current_sec.length();
        unordered_map<string, vector<string>> temp_kmerbarcodes;          

        #pragma omp parallel for schedule(guided)  shared (hetero::hetero_kmers, hetero::kmerbarcodes, current_sec, str_length)
        for (unsigned int x=0; x < (str_length-((opt::k)-1)); x++) {
            string it = current_sec.substr(x,(opt::k));
            //pos=pos+1;
            //cout << "it = " << it << endl;
            if((hetero::hetero_kmers.find(canonical_representation(it.c_str())) != hetero::hetero_kmers.end()) && (rec.comment!="")) {
            //if(bf.contains(it.c_str())) {
                //std::cout << "Present kmer " << it << " in "  << rec.id << " Barcode = "<< rec.comment <<"\n";
                #pragma omp critical
                temp_kmerbarcodes[it].push_back(rec.comment);
            }
            //A true SNP position will have multiple kmers             
        }
        
        if(temp_kmerbarcodes.size()<((opt::k)+5) && temp_kmerbarcodes.size()>((opt::k)-5)) {
            for (auto e = temp_kmerbarcodes.begin(); e != temp_kmerbarcodes.end(); e++)    {
                 for(auto &d: e->second){hetero::kmerbarcodes[e->first].push_back(d);}
            }
        //for (auto e = hetero::kmerbarcodes.begin(); e != hetero::kmerbarcodes.end(); e++)    {
        //         cout << "e->first = " << e->first << "-->"; for(auto &d: e->second){cout << "d in e->second = " << d << ",";} cout << endl;
        //    }
        temp_kmerbarcodes.clear(); 
        }

    }
    cerr << "\nFinished findcommonbarcodes for loop\n" << "hetero::kmerbarcodes size = "<< hetero::kmerbarcodes.size() << endl;   

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


Graph reduce(Graph graph) {

    /* vertex iterator */
    using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;

    //std::cout << "# of vertices: " << num_vertices(graph) << "\n";
    //std::cout << "# of edges:    " << num_edges(graph)    << "\n";

    std::set<vertex_descriptor> to_remove;

    /* iterator throught the graph  */
    for (auto self : mir(vertices(graph)))
    {      
        //cout << "self = " << self << endl; 
        //Iterate through all edges belonging to a vertex and remove edges less than the given weight
        for(auto edge : mir(out_edges(self, graph))) {
            auto weight    = boost::get(boost::edge_weight, graph, edge);
            if (weight < 10.0f) {
                //cout << "edge = " << edge << endl;;
                remove_edge(edge, graph);
            }
        }
        //cout << "degree of vertex =  " << degree(self, graph)  << "\n";
        //Remove vertex if its degree is zero 
        if (degree(self, graph) == 0)
            remove_vertex(self, graph);        
               
    }

    std::cout << "# of vertices: " << num_vertices(graph) << "\n";
    std::cout << "# of edges:    " << num_edges(graph)    << "\n";

    return graph;
}



void makegraph(){
    Graph graph;

    int vertex_id = 0;
    unordered_map<string,int> barcodeID;
    //boost::graph_traits<Graph>::vertex_descriptor vertexA = vertex(0,graph);
    //boost::graph_traits<Graph>::vertex_descriptor vertexB = vertex(0,graph);
    boost::graph_traits<Graph>::vertex_descriptor vertexA = add_vertex(Vertex {0, "DummyA" }, graph);
    boost::graph_traits<Graph>::vertex_descriptor vertexB = add_vertex(Vertex {1, "DummyB"}, graph);
    //remove_vertex(vertexA, graph);
    //remove_vertex(vertexB, graph);
    graph.clear();


    for (auto it = hetero::kmerbarcodes.begin(); it != hetero::kmerbarcodes.end(); it++) {
        // remove duplicate barcodes from this heterozygous kmer
        std::vector<std::string>& barcodes_for_kmer = it->second;
        set<string> s( barcodes_for_kmer.begin(), barcodes_for_kmer.end() );
        barcodes_for_kmer.assign( s.begin(), s.end() );

        // skip if too few barcodes
        if(barcodes_for_kmer.size() <= 1) {
            continue;
        }

        // create vertices in first pass
        for(std::size_t barcode_idx = 0; barcode_idx < barcodes_for_kmer.size(); ++barcode_idx) {
            const char* barcode_str = barcodes_for_kmer[barcode_idx].c_str();
            if (barcodeID.find(barcode_str) == barcodeID.end()) {
                vertexA = add_vertex(Vertex {vertex_id, barcode_str}, graph);
                barcodeID[barcode_str] = vertex_id;
                vertex_id++;
            }
        }

        // create edges between pairs of barcodes
        for(std::size_t barcode_idx1 = 0; barcode_idx1 < barcodes_for_kmer.size(); ++barcode_idx1) {
            const std::string& barcode_str1 = barcodes_for_kmer[barcode_idx1];
             for(std::size_t barcode_idx2 = barcode_idx1 + 1; barcode_idx2 < barcodes_for_kmer.size(); ++barcode_idx2) {
                const std::string& barcode_str2 = barcodes_for_kmer[barcode_idx2];
                //cout <<"kmer=" << it->first << ", barcode_str1 = "  << barcode_str1  << ", barcodeID[barcode_str1] ="<< barcodeID[barcode_str1] 
                //     << ", barcode_str2 = "  << barcode_str2 << ", barcodeID[barcode_str2] ="<< barcodeID[barcode_str2] << endl;
                assert(barcode_str1 != barcode_str2);
                std::pair<Edge, bool> edge_iter = (boost::edge(barcodeID[barcode_str1], barcodeID[barcode_str2], graph));

                // was edge created?
                if(edge_iter.second) {
                    int weight = get(boost::edge_weight_t(), graph, edge_iter.first);
                    //cerr << "Weight = " << weight << endl;
                    int weightToAdd = 1.0;
                    boost::put(boost::edge_weight_t(), graph, edge_iter.first, weight+weightToAdd);
                }

                if (!(edge_iter.second) 
                     && (barcodeID[barcode_str1])!=(barcodeID[barcode_str2])) {
                     //cout << "New edge\n";                     
                     //cout << "Vertex A= " << graph[barcodeID[barcode_str1]].id <<  graph[barcodeID[barcode_str1]].name << endl;
                     //cout << "Vertex B= " << graph[barcodeID[barcode_str2]].id <<  graph[barcodeID[barcode_str2]].name << endl;
                     vertexA = boost::vertex(barcodeID[barcode_str1], graph);
                     vertexB = boost::vertex(barcodeID[barcode_str2], graph);
                     add_edge(vertexA, vertexB, 1.0f, graph);
                }
             }
         }
    }

      
    std::cout << "# of vertices : " << num_vertices(graph) << "\n";
    std::cout << "# of edges:    " << num_edges(graph)    << "\n";
       
    for (auto vd : boost::make_iterator_range(vertices(graph))) {
    std::cout << "Vertex descriptor #" << vd 
         << " degree:" << degree(vd, graph)
         << " id:"     << graph[vd].id
         << " name:"     << graph[vd].name
         << "\n";
    }
    
    save(graph, "before.dot");
    numcomponents(graph);
  
    cerr << "Removing edges with weights less than 10" << endl;
    auto const h = reduce(graph);
    save(h, "after.dot");
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
        case 'k':
            arg >> opt::k;
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
/*
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

        //BloomFilter bloom(opt::m, opt::nhash, opt::k, opt::inputBloomPath.c_str());
        //cerr << bloom.getPop() << "\n";
    }
    */
    
    if (!opt::inputBloomPath.empty()) {
        findheteroreads(optind, argv);
        findcommonbarcodes(optind, argv);
        makegraph();
    }
    else
        cerr << "No input file\n";
}


/******************************************************
Code Written and Developed @ Loyola Unviersity Chicago
Jonathon Brenner, Catherine Putonti 2015
www.putonti-lab.com/software
******************************************************/

#include <vector>
#include <bitset>
#include <cstring>
#include <fstream>
#include <iostream>
#include "stdint.h"

using namespace std;

class GenomeHash{
private:
    int word_size;
    unsigned long long max_size;
    int subj_size;
    int query_size;
    int num_mm;
    int min_match_size_threshold;
    double threshold;
    typedef unsigned long long UINT64;
    static const UINT64 ONE=1;
    struct NodeU{char word[13];};//string word;};
    struct NodeS{bool num_occ=0; vector<NodeS*> MM; vector<NodeU*> ptr;};
    struct NodeO{string header;};
    bool remove_submatches();
    NodeU * htU;                // list by appearance
    NodeS * htS;                // list by hash
    vector <NodeO> results;
    vector<string> words;       // vector of sequence as words
    char * file_name;
    char * query_name;
public:
    GenomeHash(char * file_path, char * query_path, int kmer, int mm, float u_threshold);
    bool get_words_host();
    void Create();
    void Hash();
    unsigned long long word_to_decimal(string word);
    bool compare_query_multi();
    unsigned long long ptr_to_index(unsigned long long word, int pointer_index);
    unsigned long long get_match(string header, unsigned long long current_word, uint64_t i, int mm, char strand, bool &hit); // pass in query index, word; return new index
    void Get_Results();
    bool check_word(unsigned long long index);
    void get_reverse_complement();
 };

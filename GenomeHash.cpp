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
#include <cstddef>
#include "GenomeHash.h"
#include "stdint.h"
#include <algorithm>
#include <iterator>

using namespace std;

/************************************ CONSTRUCTOR ************************************
Functionality:
    Sets up the arrays
Parameters:
    file_path: name of the file for which you are screening for (intended target)
    query_path: name of the file of reads
    kmer: word size
    mm: number of mismatches in seed to consider
    u_threshold: user defined threshold (between 0 and 1)
Return:
    NULL
*************************************************************************************/
GenomeHash::GenomeHash(char * file_path, char * query_path, int kmer, int mm, float u_threshold)
{
    query_name = query_path;
    file_name = file_path;
    word_size = kmer;
    num_mm = mm;
    threshold=u_threshold;

    max_size=1;
    for(int i=0; i<word_size; i++) max_size=max_size*4;

    cout << "Reading: " << file_name << endl;
    if(get_words_host())
        Create();
    else
        cout<<"Error reading target genome"<<endl;
}

/*************************************************************************************
Functionality:
    Gets the words from the host file
Parameters:
    None
Return:
    True if correctly read
*************************************************************************************/
bool GenomeHash:: get_words_host()
{
    char * word=new char [word_size+1];
    int i,j;
    unsigned int len;
    ifstream in;
    in.open(file_name);
    if(!in.is_open()) {
        cout << "The file could not be opened. Check the location.\n";
        return false;
    }

    //gets the length
    string str_word,line,header;
    getline(in,header);       //gets header
    len=0;
    while(in.peek()!=EOF)
    {
        getline(in,str_word);
        if(str_word[0]!='>')
            len+=str_word.size();
    }
    cout << "len=" << len << endl;
    subj_size=len;
    in.clear();
    in.close();

    htU = new NodeU [len];
    unsigned long long max_size=1;
    for(int i=0; i<word_size; i++)
        max_size=max_size*4;     //calculates = 4^word_size
    htS = new NodeS [max_size];

    int ii=0;

    in.open(file_name);
    getline(in,header);       //gets header
    cout << header << endl;
    len=0;
    char c;
    while(in.peek()!=EOF)
    {
        for(i=0; i<word_size; i++)
        {
            in>>c;
            if(c<97) c+=32;
            htU[ii].word[i]=c;
        }
        ii++;

        str_word="";
        while(!(in.peek()=='>' || in.peek()==EOF))
        {
            getline(in,line);
            str_word+=line;
            for(i=0;i<(str_word.length()-word_size+1);i++)
            {
                for(j=0;j<word_size-1;j++)
                    htU[ii].word[j]=htU[ii-1].word[j+1];
                c=str_word[i];
                if(c<97) c+=32;
                htU[ii].word[word_size-1]=c;
                ii++;
            }
            //adjust str_word
            j=str_word.length()-word_size+1;
            if(!(in.peek()=='>' || in.peek()==EOF || j<(word_size-1))) str_word.erase(0,j);
        }

        if(in.peek()!=EOF)
        {
            getline(in,header);       //gets header
            cout << header << endl;
        }
    }
    subj_size=ii;
    in.clear();
    in.close();
    return true;
}

/*************************************************************************************
Functionality:
    Calls hash function and comparisons; reports to console # results
Parameters:
    None
Return:
    NULL
*************************************************************************************/
void GenomeHash::Create()
{
    Hash();
    compare_query_multi();
    cout << results.size() << " results found." << endl;
}

/*************************************************************************************
Functionality:
    Populates the htS array
Parameters:
    None
Return:
    NULL
*************************************************************************************/
void GenomeHash::Hash()
{
    unsigned long long hash_val;
    for (int i = 0; i <subj_size; i++)
    {
        hash_val=word_to_decimal(htU[i].word);
        if(! (hash_val==0 && htU[i].word[0]!='a'))
        {
            htS[hash_val].num_occ=true;
            htS[hash_val].ptr.push_back(&htU[i]);
        }
    }
}

/*************************************************************************************
Functionality:
    Hash function, converts word to integer
Parameters:
    word: nucleotide sequence
Return:
    integer value of word
*************************************************************************************/
unsigned long long GenomeHash::word_to_decimal(string word)
{
    int word_size=word.length();
    unsigned long long val = 0; // 4 bytes =  32 bits = 16 chars
    char c;
    unsigned char temp;
    int i=0;

    c=word.at(i);

    switch(c)
    {
    case 'a':
        temp=0;
        break; //shifts 00 for a
    case 't':
        temp=1;
        break; //shifts 10 for t
    case 'c':
        temp=2;
        break;//shifts 01 for c
    case 'g':
        temp=3;
        break;//shifts 11 for g
    default: return 0;
    }
    val=temp;

    for (i=1; i<word_size; i++)
    {
        val=val<<2;
        c=word.at(i);
        switch(c)
        {
        case 'a':
            temp=0;
            break; //shifts 00 for a
        case 't':
            temp=1;
            break; //shifts 10 for t
        case 'c':
            temp=2;
            break;//shifts 01 for c
        case 'g':
            temp=3;
            break;//shifts 11 for g
        default: return 0;
        }
        val+=(int)temp;
    }

    if (val>max_size)
        cout <<"Out of range: " + word <<endl;

    return val;
}

/*************************************************************************************
Functionality:
    Compares the query/read sequence(s) to the target arrays
Parameters:
    None
Return:
    True if completed correctly; false if the file cannot be opened
*************************************************************************************/
bool GenomeHash:: compare_query_multi()
{
    ifstream in;
    in.open(query_name);
    if(!in.is_open())
    {
        cout << "The read file could not be opened. Check the location.\n";
        return false;
    }
    int i,j,len;
    char * word=new char [word_size+1];

    cout << "Reading: " << query_name << endl;
    string str_word,line,header;
    getline(in,header);       //gets header

    while(in.peek()!=EOF)
    {
        len=0;
        str_word="";
        while( !(in.peek()=='>' || in.peek()==EOF))
        {
            getline(in,line);
            len+=line.size();
            str_word+=line;
        }

        min_match_size_threshold=len*threshold;


        for(i=0; i<word_size; i++)
        {
            word[i]=str_word[i];
            if(word[i]<97) word[i]+=32;                         //makes it lowercase
        }
        word[word_size]='\0';
        words.push_back(word);

        for(i=1; i<(len-word_size+1); i++)                      //read until the end of the file
        {
            //shift
            for(j=0; j<(word_size-1); j++) word[j]=word[j+1];
            word[(word_size-1)]=str_word[word_size+i-1];
            if(word[word_size-1]<97) word[word_size-1]+=32;     //makes it lowercase
            word[word_size]='\0';
            words.push_back(word);
        }

        query_size=words.size();
        bool hit=false;

        for (uint64_t ii=0; ii<query_size;ii++)                    // forward strand of query
        {
            UINT64 current_word = word_to_decimal(words[ii]);
            ii = get_match(header, current_word, ii, num_mm, '+',hit);
            if(hit) {ii=query_size;}
        }

        if(!hit)
        {
            //consider reverse complement of the query sequence
            get_reverse_complement();

            for (uint64_t ii =0; ii<query_size;ii++)                   // reverse strand of query
            {
                UINT64 current_word = word_to_decimal(words[ii]);
                ii = get_match(header, current_word, ii, num_mm, '-',hit);
                if(hit) {ii=query_size;}
            }
        }

        words.clear();

        if(in.peek()!=EOF)
        {
            getline(in,header);       //gets header
        }
    }

    in.clear();
    in.close();

    Get_Results();

    return true;
 }

/*************************************************************************************
Functionality:
    Calculates the location of the word in unsorted table of target genome
Parameters:
    word: the integer representation of the nucleotide word
    pointer_index: which particular pointer in the array of pointers
Return:
    True if completed correctly; false if the file cannot be opened
*************************************************************************************/
unsigned long long GenomeHash:: ptr_to_index(uint64_t word, int pointer_index)
{
    ptrdiff_t z;
    NodeU * base = htU;
    NodeU * s=htS[word].ptr[pointer_index];
    z=s-base;
    return (uint64_t)z;
}

/*************************************************************************************
Functionality:
    Identifies matches between target and reads meeting the threshold
Parameters:
    header: the sequence being considred
    current_word: the seed word
    i: index
    mm: number of mismatches
    strand: +/- strand being considered
    hit: boolean value passed by reference if a match is found
Return:
    index in sequence to next be considered if match found or not
*************************************************************************************/
unsigned long long GenomeHash:: get_match(string header, unsigned long long current_word, uint64_t i, int mm, char strand, bool &hit)
{ // pass in query index, word; return new index
    string match;
    double mismatches;
    double percent_mm=0;
    double match_size = 0;
    uint64_t start;
    uint64_t ii=0;
    uint64_t min_ii=query_size;
    hit=false;
    if(current_word<0 || current_word>max_size) return 0;

    if (htS[current_word].num_occ)
    {
        //loop through for each instance of query word in the subject sequence
        for(int np=0; np<htS[current_word].ptr.size(); np++)
        {
            NodeO temp;            uint64_t index = ptr_to_index(current_word,np);        // calculate index of current word in unsorted table of target genome
            match = htU[index].word;
            mismatches=0;
            percent_mm=0;
            start=i;
            ii=i;
            //ii is position within the query (words), index is position within the subject

            //jump start => match=2*word_size
            while (((ii-start)<=word_size) && (ii<(query_size-1)) && (index<(subj_size-1)))
            {
                ii++;                                           // increment query position
                index +=1;                                      // increment target position
                match+=htU[index].word[word_size-1];            // add next char to match
                if (words[ii]!=htU[index].word)
                    mismatches++;
            }
            match_size=match.size();
            percent_mm = mismatches/match_size;

            //if jump started match meets the threshold, continue to grow the word
            if(percent_mm<=threshold)
            {                while ((percent_mm<=threshold) && (ii<(query_size-1)) && (index<(subj_size-1)))
                {
                    ii++;                                       // increment query position
                    index +=1;                                  // increment target position
                    match+=htU[index].word[word_size-1];        // add next char to match
                    match_size = match.size();
                    if (words[ii]!=htU[index].word)
                    {
                        mismatches++;
                        percent_mm = mismatches/match_size;
                    }
                }

                //if the match meets the user specified threshold (for length)
                if (match_size>=min_match_size_threshold)
                {
                    temp.header=header;
                    results.push_back(temp);
                    hit=true;
                    return ii;
                }
                if(ii<min_ii) min_ii=ii;
            }
            else
                return i;
        }
    }
    else
        return i;

    return min_ii;
}

/*************************************************************************************
Functionality:
    Writes out results
Parameters:
    None
Return:
    NULL
*************************************************************************************/
void GenomeHash:: Get_Results()
{
    ofstream out_map; ofstream out_dont;
    char file[1000];
    sprintf(file,"%s.map",query_name);
    out_map.open(file);
    sprintf(file,"%s.unmapped",query_name);
    out_dont.open(file);
    ifstream in;
    in.open(query_name);
    string line;
    int i=0;
    while(in.peek()!=EOF)
    {
        getline(in,line);
        if(line[0]=='>')    //new read
        {
            if(line==results[i].header)    //mapped
            {
                out_map<<line<<endl;
                while(!(in.peek()=='>'||in.peek()==EOF))
                {
                    getline(in,line);
                    out_map<<line<<endl;
                }
                i++;
            }
            else
            {
                out_dont<<line<<endl;
                while(!(in.peek()=='>'||in.peek()==EOF))
                {
                    getline(in,line);
                    out_dont<<line<<endl;
                }
            }
        }

    }

    out_map.clear(); out_map.close();
    out_dont.clear(); out_dont.close();
    in.clear(); in.close();
}

/*************************************************************************************
Functionality:
    Checks if the word occurs within the htS array
Parameters:
    index: index in htS to be checked
Return:
    NULL
*************************************************************************************/
bool GenomeHash::check_word(unsigned long long index){
    if (htS[index].num_occ)
        return true;
    else return false;
}

/*************************************************************************************
Functionality:
    Creates the vector of words for the complementary strand of the sequence
Parameters:
    None
Return:
    NULL
*************************************************************************************/
void GenomeHash:: get_reverse_complement()
{
    vector<string> rc;
    string rev_word;
    for(int j=query_size-1; j>=0; j--)
    {
        rev_word=words[j];
        reverse(rev_word.begin(), rev_word.end());
        for(int i=0; i<rev_word.size(); i++)
        {
            switch(rev_word[i]){
            case 'a': rev_word[i]='t'; break;
            case 't': rev_word[i]='a'; break;
            case 'c': rev_word[i]='g'; break;
            case 'g': rev_word[i]='c'; break;
            default: rev_word[i]='n';
            }
        }
        rc.push_back(rev_word);
    }
    words.clear();
    words.swap(rc);
}

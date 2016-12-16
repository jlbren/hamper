/******************************************************
Code Written and Developed @ Loyola Unviersity Chicago
Jonathon Brenner, Catherine Putonti 2015
www.putonti-lab.com/software
******************************************************/

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "GenomeHash.h"

using namespace std;

/*** INPUT:
1=host file
2=query file
3=word_size
4=threshold
***/

int main(int argc, char * argv[])
{
    clock_t t1, t2;
    t1 = clock();

    char * host_file=new char [strlen(argv[1])+1];
    sprintf(host_file,"%s",argv[1]);
    char * query_file=new char [strlen(argv[2])+1];
    sprintf(query_file,"%s",argv[2]);
    int word_size=atoi(argv[3]);        //set maximum to 10. This is memory dependent. Assuming 16GB RAM. Adjust accordingly.
    if(word_size>10) word_size=10;
    if(word_size<=3) word_size=8;
    float threshold=atof(argv[4]);
    if(threshold<0.0 || threshold>1.0) threshold=0.5;
    int mm=0;

    GenomeHash * Oddish = new GenomeHash(host_file, query_file, word_size, mm, threshold);
    t2 = clock();
    float diff = t2-t1;
    float secs = diff/CLOCKS_PER_SEC;

    cout<<'\n'<<"#RUN TIME: "<<secs;

    return 0;

}

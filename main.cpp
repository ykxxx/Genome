//
//  main.cpp
//  Project 4
//
//  Created by ykx on 3/9/19.
//  Copyright © 2019 UCLA. All rights reserved.
//

#include "Trie.h"
#include "provided.h"
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

void print(vector<int> results) {
    for (int i = 0; i < results.size(); i++) {
        cout << results[i] << ", ";
    }
    cout << endl;
}

void testTrie() {
    Trie<int> t;
    t.insert("", 2);
    t.insert("ayylmoa", 4);
    t.insert("ayylmaa", 7);
    cout << "Seems legit " << endl;
    vector<int> results = t.find("ayylmao", false);
    print(results);
    cout << "found" << endl;
}

vector<Genome> testGenome()
{
    vector<Genome> vg;
    // Specify the full path and name of the gene data file on your hard drive.
    string filename = "/Users/ykx/Desktop/CS32/Project/Project4/Gee-nomics/data/test.txt";
    // Open the data file and get a ifstream object that can be used to read its
    // contents.
    ifstream strm(filename);
    if (!strm)
    {
        cout << "Cannot open " << filename << endl;
        return vg;
    }    bool success = Genome::load(strm, vg); // Load the data via the stream.
    if (success)
    {
        cout << "Loaded " << vg.size() << " genomes successfully:" << endl;
        for (int k = 0; k != vg.size(); k++) {
            cout << vg[k].name() << endl;
            string frag;
            if (vg[k].extract(0, 5, frag)) {
                cout << "Fragment: " << frag << endl;
            }
            
        }
    }
    else
        cout << "Error loading genome data" << endl;
    return vg;
} // destructor for ifstream closes the file

void testGenomeMatcher(const vector<Genome>& genomes) {
    GenomeMatcher genematcher(3);
    for (int i = 0 ; i < genomes.size(); i++) {
        genematcher.addGenome(genomes[i]);
    }
    
    std::vector<DNAMatch> temp;
    bool result = genematcher.findGenomesWithThisDNA("ATGCC", 5, true, temp);
    
    if (result) {
        cout << "result: true, temp: " << endl;
        for (int i = 0; i < temp.size(); i++) {
            cout << temp[i].genomeName << " of length " << temp[i].length
            << " at position " << temp[i].position << endl; 
        }
    }
    else {
        cout<< "result: false" << endl;
    }
    
    // test related genes
    vector<GenomeMatch> genome_match;
    genematcher.findRelatedGenomes(genomes[0], 3, true, 10, genome_match);
    
    for (int i = 0; i < genome_match.size(); i++) {
        cout << "match genome " << i << ": " << endl;
        cout << "name: " << genome_match[i].genomeName << endl;
        cout << "match percentage: " << genome_match[i].percentMatch << endl;
    }
}

int main() {
    
    vector<Genome> genomes;
    genomes = testGenome();
    testGenomeMatcher(genomes);
    
    return 0;
}

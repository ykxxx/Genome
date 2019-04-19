#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    ~GenomeMatcherImpl();
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    struct GenePosition {
        GenePosition(int gene, int pos) : m_gene(gene), m_pos(pos) {}
        int m_gene;
        int m_pos;
    };
    
    // private functions
    
    int m_minSeacherLength;
    vector<Genome>* m_geneLib;
    Trie<GenePosition>* m_geneTrie;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minSeacherLength = minSearchLength;
    m_geneLib = new vector<Genome>;
    m_geneTrie = new Trie<GenePosition>;
}

GenomeMatcherImpl::~GenomeMatcherImpl() {
    // delete m_geneLib
    for (int i = 0; i < m_geneLib->size(); i++) {
        m_geneLib->erase(m_geneLib->begin());
    }
    m_geneLib->clear();
    delete m_geneLib;
    delete m_geneTrie;
}

// O(L* N)
void GenomeMatcherImpl::addGenome(const Genome& genome)
{ // O(L * N) - minSearchLength * sequenceLegnth
    int start = 0;
    int index = (int)m_geneLib->size();
    int len = genome.length();
    while (start + m_minSeacherLength <= len) {
        string key;
        const GenePosition genePos(index, start);
        genome.extract(start, m_minSeacherLength, key);
        m_geneTrie->insert(key, genePos);
        start++;
    }
    m_geneLib->push_back(genome);
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSeacherLength;
}

// O(H*F) time, where F is the length of fragment, and H is the
//number of distinct hits across all genomes where the prefix of length minSearchLength of
//fragment (or a SNiP of the fragment, if exactMatchOnly is false) was found
bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    while (!matches.empty()) {
        matches.erase(matches.begin());
    }
    int tolerance = 0;
    if (!exactMatchOnly) {
        tolerance = 1;
    }
    if ((fragment.length() < minimumLength) || (minimumLength < m_minSeacherLength)) {
        return false;
    }
    vector<GenePosition> exist = m_geneTrie->find(fragment.substr(0, minimumSearchLength()), exactMatchOnly);
    if (exist.empty()) {
        return false;
    }
    vector<DNAMatch> temp;
    for (int i = 0; i < exist.size(); i++) {
        int match_len = 0;
        int miss = 0;
        string start;
        int gene_index = exist[i].m_gene;
        int gene_pos = exist[i].m_pos;
        Genome current = (*m_geneLib)[gene_index];
        current.extract(gene_pos, 1, start);
        if (start != fragment.substr(0, 1)) {
            continue;
        }
        gene_pos++;
        match_len++;
        for (int j = 1; j < fragment.length(); j++) {
            current.extract(gene_pos, 1, start);
            string substr = fragment.substr(j, 1);
            if ((start != substr)) {
                miss++;
                if (miss > tolerance) {
                    break;
                }
            }
            gene_pos++;
            match_len++;
        }
        DNAMatch newmatch;
        newmatch.genomeName = (*m_geneLib)[gene_index].name();
        newmatch.length = match_len;
        newmatch.position = exist[i].m_pos;
        temp.push_back(newmatch);
    }
    
    unordered_map<string, int> dna_count;
    bool meet_min_len = false;
    for (int i = 0; i < temp.size(); i++) {
        if (temp[i].length < minimumLength) {
            continue;
        }
        meet_min_len = true;
        string name = temp[i].genomeName;
        auto p = dna_count.find(name);
        if (p == dna_count.end()) {
            dna_count.insert(pair<string, int>(name, temp[i].length));
        }
        else if (temp[i].length > dna_count.at(name)) {
            dna_count.at(name) = temp[i].length;
        }
    }
    if (!meet_min_len) {
        return false;
    }
    for (int i = 0; i < temp.size(); i++) {
        string name = temp[i].genomeName;
        auto p = dna_count.find(name);
        if ((p != dna_count.end()) && (temp[i].length == dna_count.at(name))) {
            matches.push_back(temp[i]);
            dna_count.erase(p);
        }
    }
    return true;
}

// O(Q * X) time, where Q is the length in DNA bases of the query
//sequence (e.g., 3 million bases), and X is the function in the big-O of your
//findGenomesWithThisDNA() method.
bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    int num_seq = query.length() / fragmentMatchLength;
    double percent_match;
    bool has_match;
    string sequence;
    vector<DNAMatch> matches;
    unordered_map<string, int> genome_count;
    
    for (int i = 0; i < num_seq; i++) {
        has_match = false;
        query.extract(i * fragmentMatchLength, fragmentMatchLength, sequence);
        has_match = findGenomesWithThisDNA(sequence, fragmentMatchLength, exactMatchOnly, matches);
        if (has_match) {
            for (int i = 0; i < matches.size(); i++) {
                auto p = genome_count.find(matches[i].genomeName);
                if (p == genome_count.end()) {
                    genome_count.insert(pair<string, int>(matches[i].genomeName, 1));
                }
                else {
                    genome_count.at(matches[i].genomeName)++;
                }
            }
        }
    }
    
    if (genome_count.size() == 0) {
        return false;
    }
    else {
        for (auto p = genome_count.begin(); p != genome_count.end(); p++) {
            percent_match = p->second * 100 / num_seq;
            if (percent_match >= matchPercentThreshold) {
                GenomeMatch new_match;
                new_match.genomeName = p->first;
                new_match.percentMatch = percent_match;
                results.push_back(new_match);
            }
        }
    }
    return true;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}

#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
}

// O(N)
bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
    char firstchar, newChar;
    string name, line, sequence;
    bool endline = false;
    bool hasBaseline = false;
    bool emptyline = false;
    genomeSource.get(firstchar);
    if (firstchar != '>') {
        return false;
    }
    while (!genomes.empty()) {
        genomes.erase(genomes.begin());
    }
    if (!getline(genomeSource, name) || name.empty()) {
        return false;
    }
    while (genomeSource) {
        genomeSource.get(newChar);
        switch (toupper(newChar)) {
            case '>':
                if (sequence.empty() || name.empty() || emptyline) {
                    return false;
                }
                genomes.push_back(Genome(name, sequence));
                getline(genomeSource, name);
                hasBaseline = false;
                sequence = "";
                break;
            case 'A':
            case 'T':
            case 'G':
            case 'C':
            case 'N':
                if (emptyline) {
                    return false;
                }
                sequence += newChar;
                endline = false;
                hasBaseline = true;
                break;
            case '\n':
                if (endline) {
                    emptyline = true;
                }
                endline = true;
                break;
        }
    }
    if (!hasBaseline) {
        return false;
    }
    genomes.push_back(Genome(name, sequence));
    return true;
}

int GenomeImpl::length() const
{
    return m_sequence.length();
}

string GenomeImpl::name() const
{
    return m_name;
}

// O(S)
bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if (position + length > m_sequence.length()) {
        return false;
    }
    
    fragment = m_sequence.substr(position, length);
        
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}

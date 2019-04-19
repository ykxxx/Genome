#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

using namespace std;

class GenomeMatcher;

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const string& key, const ValueType& value);
    vector<ValueType> find(const string& key, bool exactMatchOnly) const;

      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct ChildNode;
    
    struct TrieNode {
        TrieNode() {
            m_values = new vector<ValueType>;
            m_children = new vector<ChildNode>;
        }
        ~TrieNode() {
            for (int i = 0; i < m_values->size(); i++) {
                m_values->erase(m_values->begin());
            }
            m_values->clear();
            delete m_values;
            delete m_children;
        }
        vector<ValueType>* m_values;
        vector<ChildNode>* m_children;
    };
    
    struct ChildNode {
        ChildNode(TrieNode* node, char key) : m_key(key), m_node(node) {}
        char m_key;
        TrieNode* m_node;
    };
    
    TrieNode* m_root;
    
    // private functions
    void deleteNodes(TrieNode* current);
    vector<ValueType> findSequence(TrieNode* current, const string& key, bool exactMatchOnly) const;
};

// O(1)
template<typename ValueType>
Trie<ValueType>::Trie() {
    m_root = new TrieNode;
}

//O(N)
template<typename ValueType>
Trie<ValueType>::~Trie() {
    deleteNodes(m_root);
}

template<typename ValueType>
void Trie<ValueType>::deleteNodes(TrieNode* current) {
    int size = current->m_children->size();
    for (int i = 0; i < size; i++) {
        deleteNodes((*current->m_children)[i].m_node);
    }
    delete current;
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value) {
    // O(LC) / O(L)
    bool found = false;
    TrieNode* current = m_root;
    for (int i = 0; i < key.length(); i++) {
        found = false;
        for (int j = 0; j < current->m_children->size(); j++) {
            // if the key has been found, continue with the key
            if ((!current->m_children->empty()) && ((*current->m_children)[j].m_key == key[i])) {
                current = (*current->m_children)[j].m_node;
                found = true;
                break;
            }
        }
        // if the no node with key exists, create new
        if (!found) {
            TrieNode* newOne = new TrieNode;
            current->m_children->insert(current->m_children->end(), ChildNode(newOne, key[i]));
            current = newOne;
        }
    }
    current->m_values->push_back(value);
}

// true O(LC+V)
// false O(L2 C2+V)
// where L is the
//length of the searched-for key, C is the average number of children per node in your trie,
//and V is the size of the returned vector.
template<typename ValueType>
vector<ValueType> Trie<ValueType>::find(const string& key, bool exactMatchOnly) const {
    vector<ValueType> results;
    vector<ValueType> r;
    TrieNode* current = m_root;
    for (int i = 0; i < current->m_children->size(); i++) {
        //while (!(*current->m_children).empty()) {
            if ((*current->m_children)[i].m_key == key[0]) {
                r = findSequence((*current->m_children)[i].m_node, key.substr(1), exactMatchOnly);
                //r = findSequence(current, key, 0, exactMatchOnly);
                results.insert(results.end(), r.begin(), r.end());
            }
            //current = (*current->m_children)[i].m_node;
        //}
    }
    return results;
}

template<typename ValueType>
vector<ValueType> Trie<ValueType>::findSequence(TrieNode* node, const string& key, bool exactMatchOnly) const {
    bool found = false;
    vector<ValueType> result;
    if (key.length() == 0) {
        return *(node->m_values);
    }
    int size = node->m_children->size();
    for (int i = 0; i < size; i++) {
        ChildNode current = (*node->m_children)[i];
        char c_key = current.m_key;
        char n_key = key[0];
            if ((c_key == n_key)){
                found = true;
                result = findSequence(current.m_node, key.substr(1), exactMatchOnly);
            }
        if (!result.empty()) {
            break;
        }
    }
    if (!found) {
        if (exactMatchOnly) {
            return result;
        }
            exactMatchOnly = true;
            for (int i = 0; i < size; i++) {
                ChildNode current = (*node->m_children)[i];
                result = findSequence(current.m_node, key.substr(1), exactMatchOnly);
        }
    }
    return result;
}

//O(N)
template<typename ValueType>
void Trie<ValueType>::reset() {
    deleteNodes(m_root);
    TrieNode* new_root = new TrieNode;
    m_root = new_root;
}

#endif // TRIE_INCLUDED

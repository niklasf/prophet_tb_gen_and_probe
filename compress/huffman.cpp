#include "huffman.h"
#include <vector>
#include <queue>

struct NodeCmp {
    bool operator()(const Node* a, const Node* b) const {
        return a->freq > b->freq;
    }
};
Node* buildTree(uint64_t* freq, uint16_t size) {
    std::priority_queue<Node*, std::vector<Node*>, NodeCmp> pq;

    for (int i=0; i < size; i++)
        if (freq[i] > 0)
            pq.push(new Node(freq[i], i));

    if (pq.empty())
        return new Node(1, 0); // dummy tree if input empty

    while (pq.size() > 1) {
        Node* a = pq.top(); pq.pop();
        Node* b = pq.top(); pq.pop();
        pq.push(new Node(a->freq + b->freq, DUMMY_CODE, a, b));
    }
    return pq.top();
}
void extractLengths(Node* n, int depth, uint64_t* lengths) {
    if (!n) return;
    if (n->symbol != DUMMY_CODE) {
        lengths[n->symbol] = depth;
        return;
    }
    extractLengths(n->left,  depth + 1, lengths);
    extractLengths(n->right, depth + 1, lengths);
}
#ifndef HUFF_INCLUDED
#define HUFF_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include "repair.h"

struct Node {
    uint64_t freq;
    CODE symbol;      // >=0 for leaves, -1 for internal
    Node* left;
    Node* right;

    Node(uint32_t f, int s, Node* L=nullptr, Node* R=nullptr)
        : freq(f), symbol(s), left(L), right(R) {}
    ~Node() {
        delete left;
        delete right;
    }
};

Node* buildTree(uint64_t* freq, uint16_t size);
void extractLengths(Node* n, int depth, uint64_t* lengths);

#endif

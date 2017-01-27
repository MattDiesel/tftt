
#ifndef TFTT_TREE_H
#define TFTT_TREE_H

#include <set>

#include "tftt.h"
#include "treegroup.h"


namespace tftt {


struct Tree {
	TreeGroup* root;
	size_t ccells;

	double size[DIM];

    ~Tree();

	// For thread
	// TreeGroup* first;
	// TreeGroup* last;

	//! The set of ghost groups on this processor
	// std::set<TreeGroup*> ghosts;

};


//! The tree instance
extern Tree gtree;



struct {
    class leaf_iterator {
        CellRef cr;

        void next() {
            // if (!cr.isValid()) {
            //     return;
            // }

            if (cr.index+1 >= 2*DIM) {
                if (cr.group == gtree.root) {
                    cr = CellRef();
                    return;
                }
                cr = cr.parent();
                next();
            }
            else {
                cr.index++;
            }

            while (cr.hasChildren()) {
                cr = cr.child(0);
            }
        }
    public:
        leaf_iterator(CellRef c) : cr(c) {
            while (cr.hasChildren()) {
                cr = cr.child(0);
            }
        }

        leaf_iterator operator++() { leaf_iterator i = *this; next(); return i; }
        leaf_iterator operator++(int junk) { next(); return *this; }
        CellRef& operator*() { return cr; }
        CellRef* operator->() { return &cr; }
        bool operator==(const leaf_iterator& rhs) { return cr == rhs.cr; }
        bool operator!=(const leaf_iterator& rhs) { return !(cr == rhs.cr); }
    };

    leaf_iterator begin() { return leaf_iterator(CellRef(gtree.root, 0)); }
    leaf_iterator end() { return leaf_iterator(CellRef()); }
} leaves;


} // namespace tftt


#endif

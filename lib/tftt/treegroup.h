
#ifndef TFTT_TREEGROUP_H
#define TFTT_TREEGROUP_H


#include <array>

#include "tftt.h"
#include "treeid.h"
#include "cellref.h"
#include "treecell.h"


namespace tftt {


//! A group of adjacent sibling cells. (A quad in 2d, oct in 3d)
//! This reduces the memory requirement by not having to link all cells to their neighbours.
//!
//! TreeGroups should never be directly edited by the user.
struct TreeGroup {

	// Base Tree
	TreeGroup();
	TreeGroup(CellRef parent);
	~TreeGroup();

	//! The identifier of the group, equal to the identifier of the first cell.
	//! The id of other cells in the group can be found by adding their index to this id.
	ident_t id;

	//! The origin (top left corner in 2d) of the group.
	double origin[DIM];

	//! The cells in this group
	TreeCell cells[1<<DIM];

	// FTT

	//! Pointers to the neighbouring groups
	std::array<CellRef,DIM*2> neighbours;

	//! Pointer to the parent group.
	CellRef parent;

	// For Thread

	//! The next group in the space filling curve.
	// TreeGroup* next;

	//! The previous group in the space filling curve
	// TreeGroup* prev;

	// For parallel

	//! The index of the processor this group is active on.
	// int rank;

	bool isCopyBoundary() const;
	bool isReflectBoundary() const;
	bool isBoundary() const;
};


extern TreeGroup CopyNeighbour;
extern TreeGroup ReflectNeighbour;


} // namespace tftt


#endif

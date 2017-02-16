
#ifndef TFTT_HILBERT_H
#define TFTT_HILBERT_H


namespace tftt {


//! Returns the orientation of the hilbert child
int hilbOrient(int Dp, int ch);

//! Returns the cell index of hilbert child hch
int hilbChild(int Dp, int hch);

//! Returns the hilbert index of a cell
int hilbInvChild(int Dp, int ch);

//! Checks if a child is the last in its parent
bool hilbIsLast(int Dp, int ch);

//! Checks if a child is the first in its parent
bool hilbIsFirst(int Dp, int ch);


} // namespace TFTT


#endif

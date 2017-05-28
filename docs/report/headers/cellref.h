class CellRef {
private:
    TreeCell* _cell;
public:
    // Basic Orthotree functions
    data_t& data();                    // $\mathcal{C}\acc\textproc{data}$
    CellRef parent() const;            // $\mathcal{C}\acc\textproc{parent}$
    TreeGroup* children() const;       // $\mathcal{C}\acc\vect{\mathcal{K}}$
    CellRef child(int c) const;        // $\mathcal{C}\acc\mathcal{K}_c$

    // Core FTT functionality
    CellRef neighbour(int n) const;    // $\mathcal{C}\acc\mathcal{N}_n$

    // TFTT functionality
    int orientation() const;           // $\mathcal{C}\acc\textproc{orient}$
    CellRef next() const;              // $\mathcal{C}\acc\textproc{next}$
    CellRef prev() const;              // $\mathcal{C}\acc\textproc{prev}$
    CellRef firstChild() const;         // $\mathcal{C}\acc\textproc{firstChild}$
    CellRef lastChild() const;         // $\mathcal{C}\acc\textproc{lastChild}$

    // Parallel 
    node_t& rank();                    // $\mathcal{C}\acc\textproc{rank}$
    ident_t id() const;                // $\mathcal{C}\acc\textproc{ident}$

    // Poisson Neighbourhoods
    int poisNeighbours() const;        // $\abs{\mathcal{C}\acc\vect{\mathcal{P}}}$
    double poisCoef(int p) const;      // $\mathcal{C}\acc\beta_p$
    CellRef poisNeighb(int p) const;   // $\mathcal{C}\acc\mathcal{P}_p$
    double poisCentralCoef() const;    // $\mathcal{C}\acc\alpha$
};

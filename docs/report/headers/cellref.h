class CellRef {
private:
    TreeCell* _cell;
public:
    // Basic Orthotree functions
    data_t& data();              // $\mathcal{C}\acc\textproc{data}$
    CellRef parent();            // $\mathcal{C}\acc\textproc{parent}$
    TreeGroup* children();       // $\mathcal{C}\acc\vect{\mathcal{K}}$
    CellRef child(int c);        // $\mathcal{C}\acc\mathcal{K}_c$

    // Core FTT functionality
    CellRef neighbour(int n);    // $\mathcal{C}\acc\mathcal{N}_n$

    // TFTT functionality
    int& orientation();          // $\mathcal{C}\acc\textproc{orient}$
    CellRef& next();             // $\mathcal{C}\acc\textproc{next}$
    CellRef& prev();             // $\mathcal{C}\acc\textproc{prev}$
    CellRef firstChild();         // $\mathcal{C}\acc\textproc{firstChild}$
    CellRef lastChild();         // $\mathcal{C}\acc\textproc{lastChild}$

    // Parallel 
    node_t& rank();              // $\mathcal{C}\acc\textproc{rank}$
    ident_t id();                // $\mathcal{C}\acc\textproc{ident}$

    // Poisson Neighbourhoods
    int poisNeighbours();        // $\abs{\mathcal{C}\acc\vect{\mathcal{P}}}$
    double& poisCoef(int p);     // $\mathcal{C}\acc\beta_p$
    CellRef& poisNeighb(int p);  // $\mathcal{C}\acc\mathcal{P}_p$
    double& poisCentralCoef();   // $\mathcal{C}\acc\alpha$
};

/*
 *
 *  Created by David Bryant on 8/03/10.
 *
 */



#include <unistd.h>
#include <ctime>
#include <list>
#include "PhyLib/phylib.h"

using namespace Phylib;



/**
 Prints help and usage information
 **/
void printUsage(ostream& os, const string& madName) {
    os<<madName<<"\n\nMinimum Ancestral Deviation root finding\n";
    os<<"Usage:\n\nTo compute the MadRoot use:\n\t"<<madName<<" <treefile>\nTo run the simulation from the paper, use\n\t"<<madName<<" -SIM\n\n\n";
    os<<"The rooted trees are then output to stdout"<<endl;
    os<<endl;
}

/**
 Handles input.
 
 Returns false if there is an error.
 If the -SIM option is selected, then taxa and tree will be empty/null.
 
 //TODO Add verbose option.
 **/

bool handleInput(int argc,
                     char* argv[],
                     bool& runSimulation,
                     ifstream& is) {
    
   
    
    string madName = string(argv[0]);
    if (argc!=2) {
        printUsage(cerr,madName);
        return false;
    }
    
    string inputfile = string(argv[1]);
    
    if (inputfile=="-SIM") {
        runSimulation = true;
        return true;
    }
    runSimulation = false;
    
    
    is.open(inputfile.c_str());
    
    if (!is) {
        cerr<<"Error reading in file "<<inputfile<<"\n\n";
        printUsage(cerr,madName);
        return false;
    }
    
    
    
    
    
    
    return true;
}

/**
 tidyTree
 
 Converts tree to a binary tree (inserting zero length branches) which is rooted at the leaf
 which initially was the leftmost leaf.
 
 Any degree 2 nodes are removed.
 **/
void tidyTree(phylo<basic_newick>& tree) {
    typedef typename phylo<basic_newick>::iterator ITERATOR;
    
    //1 Reroot the tree
    reRoot(tree,tree.leftmost_leaf());
    
    
    //    cout<<"Rerooted tree"<<endl;
    //    print_newick<basic_newick>(cout, tree);
    //    cout<<endl;
    
    
    
    //2 Contract any degree two nodes except the root.
    //NOTE: DOESN'T WORK WITH INTERNAL LABELS.
    ITERATOR p = tree.leftmost_leaf();
    while (!p.root()) {
        if (!p.root() && !p.left().null() && p.left().right().null() ) {
            //p has a single child. Suppress it.
            ITERATOR nextp = p.next_post();
            p.left()->length += p->length;
            tree.graft_sibling(p,tree,p.left());
            tree.erase(p);
            p = nextp;
        }
        else
            p = p.next_post();
    }
    
    //    cout<<"Contracted tree"<<endl;
    //    print_newick<basic_newick>(cout, tree);
    //    cout<<endl;
    
    
    //3 Iterate through all non-root nodes. If it has more than 2 children, insert new nodes to make a
    // 'local' caterpillar tree with zero length branches, making the tree binary.
    for(ITERATOR p = tree.leftmost_leaf(); !p.root(); p = p.next_post()) {
        if (!p.left().null()) {
            //We know that p has at least two children, from step 2.
            ITERATOR q = p.left().right();
            while (!q.right().null()) {
                ITERATOR newNode = tree.insert_child(p); //Create new leftmost child of p
                newNode->length = 0.0;
                tree.graft_child(newNode,tree,q.right()); //Move the third child to become child of the newNode.
                tree.graft_child(newNode,tree,newNode.right()); //Move the old leftmost node to below newNode
            }
        }
    }
    
    
}

class nodeData :public basic_newick  {
public:
    nodeData() : basic_newick() {
        dist = a_uv = b_uv = 0.0;
    }
    nodeData(const basic_newick& node) : basic_newick(node) {
        dist = a_uv = b_uv = 0.0;
    }
    double dist;
    double a_uv;
    double b_uv;
};

/**
 Comnpute a_uv and b_uv for all external edges uv except for the edge incident with the root.
 
 We use the formulas, which are valid when u is a leaf
 a_{uv}  = \sum_{x \neq u} \frac{4}{(d_{ux})^2}
 b_{uv}  =  \sum_{x \neq u} \left( \frac{-4}{d_{ux}} \right)
 
 To compute these sums, we iterate over all pairs of leaves x,y, compute d_xy, and then
 update the sums for u=x and u=y.  We do this by looping over internal nodes z and then through
 pairs x,y such that lca(x,y)=z, as this enables easy calculation of d_{xy}.
 **/

void computeExternalCoefficients(phylo<nodeData>& tree) {
    typedef typename phylo<nodeData>::iterator ITERATOR;
    
    //Zero the coefficients (for safety's sake)
    for (ITERATOR z = tree.leftmost_leaf();z!=tree.root();z=z.next_post()) {
        z->a_uv = z->b_uv = 0.0;
    }
    
    //Loop through
    for (ITERATOR z = tree.leftmost_leaf();z!=tree.root();z=z.next_post()) {
        if (z.left().null())
            continue; //Skip leaves
        
        z->dist = 0.0;  //dist is the distance to z.
        //Compute the distances from z to each node in the subtree
        for(ITERATOR x = z.next_pre(z);!x.null();x = x.next_pre(z)) {
            x->dist = x.par()->dist + x->length;
        }
        
        //Now loop through leaves in the two subtrees branching off at z
        ITERATOR z1 = z.left();
        ITERATOR z2 = z.left().right();
        
        for(ITERATOR x = z1;!x.null();x = x.next_pre(z1)) {
            if (!x.left().null()) //skip internal nodes
                continue;
            for(ITERATOR y = z2;!y.null();y=y.next_pre(z2)) {
                if (y.left().null()) {
                    double dxy = x->dist + y->dist;
                    if (dxy > 0) {
                        x->a_uv += 4.0/(dxy*dxy);
                        y->a_uv += 4.0/(dxy*dxy);
                        x->b_uv += -4.0/dxy;
                        y->b_uv += -4.0/dxy;
                    }
                }
            }
            
        }
        
    }
    
    //Add distances from root.
    tree.root()->dist = 0.0;
    for(ITERATOR x = tree.root().next_pre();!x.null();x = x.next_pre()) {
        x->dist = x.par()->dist + x->length;
        if (x.left().null()) {
            double dxy = x->dist;
            if (dxy > 0) {
                x->a_uv += 4.0/(dxy*dxy);
                x->b_uv += -4.0/dxy;
            }
        }
    }
    
    
}

/**
 Compute a_uv and b_uv for all the internal nodes u, v = u.par().
 
 **/

void computeInternalCoefficients(phylo<nodeData>& tree) {
    typedef typename phylo<nodeData>::iterator ITERATOR;
    
    for( ITERATOR u = tree.leftmost_leaf();u!=tree.root();u=u.next_post()) {
        if (!u.left().null()) {
            //Internal node.
            ITERATOR u1 = u.left();
            ITERATOR u2 = u.left().right();
            
            //Compute a_uv using
            // a_uv = a_{u_1 u} + a_{u_2 u}   - \sum_{x \in U_1} \sum_{y \in U_2} \frac{8}{(d_{xy})^2}
            u->dist = 0.0;
            for (ITERATOR x = u.next_pre(u);!x.null();x=x.next_pre(u))
                x->dist = x.par()->dist+x->length;   //x->dist = d_xu.
            
            u->a_uv = u1->a_uv + u2->a_uv;
            for(ITERATOR x = u1;!x.null();x = x.next_pre(u1)) {
                if (!x.left().null())
                    continue;
                for(ITERATOR y = u2;!y.null();y = y.next_pre(u2)) {
                    if (y.left().null()) {
                        double dxy = x->dist + y->dist;
                        if (dxy > 0)
                            u->a_uv += -8.0 / (dxy * dxy);
                    }
                }
            }
            
            //compute b_uv
            u->b_uv = u1->b_uv + u2->b_uv + 2.0*(u1->length)*(u1->a_uv) + 2.0*(u2->length)*(u2->a_uv);
        }
    }
}


/**
 Suppose that v is the root. The recursion implemented in computeInternalCoefficients computes
 a_uv and b_uv.
 
 As a check, we here compute them directly from the definition. Normally this function would not be called.
 **/

void recomputeTopEdge(phylo<nodeData>& tree, double& a_uv, double& b_uv) {
    /*For this edge (V = \{v\} = root), we have
     a_uv = \sum_x 4/(d_{xv})^2
     b_uv = \sum_x (-4/d_{xv} + 8 d_{ux} / (d_{xv})^2)
     = \sum_x (-4/d_{xv} + 8 (d_xv - d_uv) / (d_{xv})^2)
     = \sum_x (4/d_{xv} - 8 d_uv) / (d_{xv})^2)
     */
    typedef typename phylo<nodeData>::iterator ITERATOR;
    
    ITERATOR u = tree.root().left();
    ITERATOR v = tree.root();
    double duv = u->length;
    
    v->dist = 0.0;   //  d_xv
    a_uv = 0.0;
    b_uv = 0.0;
    for (ITERATOR x = tree.root().next_pre();!x.null();x = x.next_pre()) {
        x->dist = x.par()->dist + x->length;
        if (x.left().null()) {
            double dxv = x->dist;
            a_uv += 4.0/(dxv*dxv);
            b_uv += 4.0/dxv - 8.0*duv/(dxv*dxv);
        }
    }
    
}
void contractZeroLength(phylo<nodeData>& tree) {
    typedef typename phylo<nodeData>::iterator ITERATOR;
    
    for(ITERATOR p = tree.leftmost_leaf();!p.null();p=p.next_post()) {
        //Contract any children of this node which have zero length branches.
        ITERATOR q = p.left();
        while(!q.null()) {
            if (q->length == 0.0) {
                ITERATOR nextq = q.right();
                tree.contract(q);
                q = nextq;
            } else
                q=q.right();
        }
        
    }
}

void findOptimal(phylo<nodeData>& tree, phylo<nodeData>::iterator& u, double& t) {
    typedef typename phylo<nodeData>::iterator ITERATOR;
    
    for(u = tree.leftmost_leaf();u!=tree.root();u = u.next_post()) {
        t = -u->b_uv/(2.0*(u->a_uv));
        if (0 < t && t < u->length) {
            //Optimum along internal branch.
            return;
        }
        
        if (t <= 0 && !u.left().null()) {
            bool stillOptimal = true;
            ITERATOR uk = u.left();
            while (!uk.null() && stillOptimal) {
                double tk = -uk->b_uv/(2.0*(uk->a_uv));
                if (tk < uk->length)
                    stillOptimal = false;
                uk=uk.right();
            }
            if (stillOptimal) {
                t=0;
                return;
            }
        }
    }
    cerr<<"Unable to find optimum"<<endl;
    exit(-2);
}



/**
 Determines the MAD root. If the location is along an edge, this edge is subdivided.
 The tree is rooted at the MAD root.
 
 Assume that tree is already rooted at a leaf (note: we might want to do this within the
 function in future).
 Assumes that leaves (and root) are labelled uniquely by 0,...,ntax-1.
 **/
void MADroot(phylo<basic_newick>& newickTree) {
    
    phylo<nodeData> tree;
    tidyTree(newickTree);
    
    //    cerr<<"Tidied tree is "<<endl;
    //    print_newick<basic_newick>(cerr, newickTree);
    //    cerr<<endl;
    
    copy(newickTree, tree); //Copy the generic structure over to the one we want.
    
    
    
    
    
    //    cerr<<"Copied tree is "<<endl;
    //    print_newick<nodeData>(cerr, tree);
    //    cerr<<endl;
    
    
    
    computeExternalCoefficients(tree);
    computeInternalCoefficients(tree);
    
    //Contract zero length branches
    contractZeroLength(tree);
    
    
    //-------
    //SANITY CHECK
    //    double a_uv,b_uv;
    //    recomputeTopEdge(tree,a_uv, b_uv);
    //    cerr<<"Top edge a_uv is "<<tree.root().left()->a_uv<<" and should be "<< a_uv<<endl;
    //    cerr<<"Top edge b_uv is "<<tree.root().left()->b_uv<<" and should be "<< b_uv<<endl;
    //    cerr<<endl;
    //--------
    
    phylo<nodeData>::iterator u;
    double t;
    findOptimal(tree,u,t);
    
    if (t==0) {
        reRoot(tree,u);
    } else {
        //Subdivide node to make root
        double duv = u->length;
        phylo<nodeData>::iterator newNode = tree.insert_child(u.par());
        newNode->length = duv - t;
        tree.graft_child(newNode, tree, u);
        u->length = t;
        reRoot(tree,newNode);
    }
    newickTree.clear();
    copy(tree, newickTree); //Copy the generic structure over to the one we want.
    tree.clear();
}

//This is where we specify what is stored at each node.

void runSimulation(ostream& os) {
    
    int nreps = 10;
    os<<"MadROOT\n\nSIMULATION OUTPUT\n";
    os<<"ntax\treplicate\tdistribution\ttime"<<endl;
    std::list<unsigned int> listOfN( { 100, 160, 250, 400, 640, 1000, 1600, 2500, 4000, 6400, 10000, 16000, 25000, 40000, 64000, 100000 } );
    for (int ntax : listOfN) {
        for (int treeDistr=0;treeDistr<2;treeDistr++) {
            for(int i=0;i<nreps;i++) {
                phylo<basic_newick> tree;
                if (treeDistr==0)
                    coalescent(tree,ntax,1.0);
                else
                    sim_uniform_phylo(tree,ntax);
                
                //Unit branch lengths
                for(phylo<basic_newick>::iterator p = tree.leftmost_leaf();p!=tree.root();p=p.next_post())
                    p->length = 1.0;
                
                clock_t begin = clock();
                MADroot(tree);
                clock_t end = clock();
                double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                tree.clear();
                
                os<<ntax<<"\t"<<(i+1)<<"\t";
                if (treeDistr==0)
                    os<<"COALESCENT\t";
                else
                    os<<"UNIFORM\t";
                os<<elapsed_secs<<endl;
            }
            
        }
    }
    
    
}




int main(int argc, char* argv[]) {
    
    vector<string> taxa;
    phylo<basic_newick> newickTree;
    bool runningSimulation = false;
    ifstream is;
    
    bool OK = handleInput(argc, argv, runningSimulation,is);
    if (!OK)
        exit(1);
    
    if (runningSimulation) {
        cerr<<"RUNNING SIMULATION"<<endl;
        runSimulation(cout);
        exit(0);
    }
    
    
    //Read in trees one at a time and output rooted versions.
    
    bool done = false;
    while(!done) {
        try {
         read_newick(is, newickTree, taxa, 0.0);
        } catch(...) {
            done = true;
        }
         if (taxa.empty()) {
             done = true;
         }
        else {
            MADroot(newickTree);
            print_newick<basic_newick>(cout, newickTree,taxa,true,true);
            cout<<";"<<endl;
        }
        taxa.clear();
        newickTree.clear();
        
    }
    
    is.close();
    
    return(0);
}

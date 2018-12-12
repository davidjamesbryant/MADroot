library("inline")
library("Rcpp")

settings <- getPlugin("Rcpp")
settings$env$PKG_CXXFLAGS <- paste("-I",getwd(),sep="")

madRootCppCode<-'


#include <unistd.h>
#include <ctime>
#include <list>
#include "Phylib/phylib.h"

using namespace Phylib;



/**
tidyTree

Converts tree to a binary tree (inserting zero length branches) which is rooted at the leaf
which initially was the leftmost leaf.

Any degree 2 nodes are removed.
**/
void tidyTree(phylo<basic_newick>& tree) {

/* 1 Reroot the tree */
reRoot(tree,tree.leftmost_leaf());



/* 2 Contract any degree two nodes except the root. */
/* NOTE: DOESNT WORK WITH INTERNAL LABELS. */
phylo<basic_newick>::iterator  p = tree.leftmost_leaf();
while (!p.root()) {
if (!p.root() && !p.left().null() && p.left().right().null() ) {
//p has a single child. Suppress it.
phylo<basic_newick>::iterator  nextp = p.next_post();
p.left()->length += p->length;
tree.graft_sibling(p,tree,p.left());
tree.erase(p);
p = nextp;
}
else
p = p.next_post();
}


//3 Iterate through all non-root nodes. If it has more than 2 children, insert new nodes to make a
// local caterpillar tree with zero length branches, making the tree binary.
for(phylo<basic_newick>::iterator  p = tree.leftmost_leaf(); !p.root(); p = p.next_post()) {
if (!p.left().null()) {
//We know that p has at least two children, from step 2.
phylo<basic_newick>::iterator  q = p.left().right();
while (!q.right().null()) {
phylo<basic_newick>::iterator  newNode = tree.insert_child(p); //Create new leftmost child of p
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



void computeExternalCoefficients(phylo<nodeData>& tree) {

//Zero the coefficients (for safetys sake)
for (phylo<nodeData>::iterator z = tree.leftmost_leaf();z!=tree.root();z=z.next_post()) {
z->a_uv = z->b_uv = 0.0;
}

//Loop through
for (phylo<nodeData>::iterator z = tree.leftmost_leaf();z!=tree.root();z=z.next_post()) {
if (z.left().null())
continue; //Skip leaves

z->dist = 0.0;  //dist is the distance to z.
//Compute the distances from z to each node in the subtree
for(phylo<nodeData>::iterator x = z.next_pre(z);!x.null();x = x.next_pre(z)) {
x->dist = x.par()->dist + x->length;
}

//Now loop through leaves in the two subtrees branching off at z
phylo<nodeData>::iterator z1 = z.left();
phylo<nodeData>::iterator z2 = z.left().right();

for(phylo<nodeData>::iterator x = z1;!x.null();x = x.next_pre(z1)) {
if (!x.left().null()) //skip internal nodes
continue;
for(phylo<nodeData>::iterator y = z2;!y.null();y=y.next_pre(z2)) {
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
for(phylo<nodeData>::iterator x = tree.root().next_pre();!x.null();x = x.next_pre()) {
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

for( phylo<nodeData>::iterator u = tree.leftmost_leaf();u!=tree.root();u=u.next_post()) {
if (!u.left().null()) {
//Internal node.
phylo<nodeData>::iterator u1 = u.left();
phylo<nodeData>::iterator u2 = u.left().right();

//Compute a_uv
u->dist = 0.0;
for (phylo<nodeData>::iterator x = u.next_pre(u);!x.null();x=x.next_pre(u))
x->dist = x.par()->dist+x->length;   //x->dist = d_xu.

u->a_uv = u1->a_uv + u2->a_uv;
for(phylo<nodeData>::iterator x = u1;!x.null();x = x.next_pre(u1)) {
if (!x.left().null())
continue;
for(phylo<nodeData>::iterator y = u2;!y.null();y = y.next_pre(u2)) {
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


phylo<nodeData>::iterator u = tree.root().left();
phylo<nodeData>::iterator v = tree.root();
double duv = u->length;

v->dist = 0.0;   //  d_xv
a_uv = 0.0;
b_uv = 0.0;
for (phylo<nodeData>::iterator x = tree.root().next_pre();!x.null();x = x.next_pre()) {
x->dist = x.par()->dist + x->length;
if (x.left().null()) {
double dxv = x->dist;
a_uv += 4.0/(dxv*dxv);
b_uv += 4.0/dxv - 8.0*duv/(dxv*dxv);
}
}

}
void contractZeroLength(phylo<nodeData>& tree) {

for(phylo<nodeData>::iterator p = tree.leftmost_leaf();!p.null();p=p.next_post()) {
//Contract any children of this node which have zero length branches.
phylo<nodeData>::iterator q = p.left();
while(!q.null()) {
if (q->length == 0.0) {
phylo<nodeData>::iterator nextq = q.right();
tree.contract(q);
q = nextq;
} else
q=q.right();
}

}
}

void findOptimal(phylo<nodeData>& tree, phylo<nodeData>::iterator& u, double& t) {

for(u = tree.leftmost_leaf();u!=tree.root();u = u.next_post()) {
t = -u->b_uv/(2.0*(u->a_uv));
if (0 < t && t < u->length) {
//Optimum along internal branch.
return;
}

if (t <= 0 && !u.left().null()) {
bool stillOptimal = true;
phylo<nodeData>::iterator uk = u.left();
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
copy(newickTree, tree); //Copy the generic structure over to the one we want.

computeExternalCoefficients(tree);
computeInternalCoefficients(tree);

//Contract zero length branches
contractZeroLength(tree);

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

bool isTreeTrivial(phylo<basic_newick>& tree) {
if (tree.root().null())
return true;
int leafCount = 0;
for(phylo<basic_newick>::iterator p = tree.root();!p.null();p=p.next_pre())
if (p.left().null())
leafCount++;
return (leafCount<=1);
}


vector<std::string> applyMADroot(std::string input) {

vector<string> taxa;
phylo<basic_newick> newickTree;
stringstream is(input);
stringstream os;
vector<std::string> outputs;


//Read in trees one at a time and output rooted versions.

bool done = false;
while(!done) {
try {
read_newick(is, newickTree, taxa, 0.0);
} catch(...) {
done = true;
}
if (taxa.empty()||isTreeTrivial(newickTree)) {
done = true;
}
else {
MADroot(newickTree);
print_newick<basic_newick>(os, newickTree,taxa,true,true);
os<<";";
outputs.push_back(os.str());
os.str("");

}
taxa.clear();
newickTree.clear();

}
return(outputs);
}

'



##nowusethesnippedaboveaswellasoneargumentconversion 
## in as well as out to provide Fibonacci numbers via C++ 
madRoot<-cxxfunction(signature(xc="character"),
          plugin="Rcpp", 
          incl=madRootCppCode, 
          settings = settings,
          body='
    std::string inputs = Rcpp::as<std::string>(xc);
    vector<std::string> outputs = applyMADroot(inputs);
    return Rcpp::wrap(outputs);
')

//DONT FORGET to merge after split
//Group ops are usable only for implicit treap, primary key - for not implicit.
//Simple made persistent by creating new verts instead of modifying existing
typedef int pkey_t; typedef int val_t;//pkey_t = key type, val_t = data key
typedef int rng_t; typedef int pval_t; //rng_t = rand key, pval_t = push value
const static val_t neutral = 0; const static pval_t neutralPush = 0;
struct Node {
    rng_t sk; int sz; //sk - rand key, sz - subtree sz
    Node *left, *right; //left,right - pointers to childs or null
    val_t val, fVal;//val - value in vert, fVal - f on subtree
    pval_t def; pkey_t pk; //def - deferred change, pk - search tree key
    Node(pkey_t key, val_t val) {
        pk = key; this->val = fVal = val; sk = rand(); sz = 1; left = right = nullptr; def = neutralPush;
    }
};
int size(Node* root) {return root ? root->sz : 0;}
val_t getF(Node* root) {return root ? root->fVal : neutral;};
val_t rqOp(val_t lhs, val_t rhs) {return lhs+rhs;}//OPERATION FOR QUERIES IN TREAP
void recalc(Node* root) {
    root->sz = 1 + size(root->left) + size(root->right);
    root->fVal = rqOp(getF(root->left), rqOp(root->val, getF(root->right)));
}
void push(Node* root) {
    //push delayed modifications to children (if exists) + maybe swap them (example: tree with range rev)
    if (root->left) root->left->def = combineUpdates(root->left->def, root->def);
    if (root->right) root->right->def = combineUpdates(root->right->def, root->def);
    root->def = neutralPush;
}
//split by pk "<" to the left, ">="(not "<") to the right
pair<Node*, Node*> split(Node* root, pkey_t key) {
    if (!root) return {nullptr, nullptr};
    push(root); auto rootVal = root->pk;
    if (rootVal < key) {//criteria
        auto p = split(root->right, key); root->right = p.first;
        recalc(root); return {root, p.second};
    } else {
        auto p = split(root->left, key); root->left = p.second;
        recalc(root); return {p.first, root};
    }
}
pair<Node*, Node*> splitKLeft(Node* root, int k) {//cut k items from left (most usable in implicit treap)
    if (!root) return {nullptr, nullptr};
    push(root); auto leftSz = size(root->left);
    if (leftSz < k) {
        auto p = splitKLeft(root->right, k-leftSz-1); root->right = p.first;
        recalc(root); return {root, p.second};
    } else {
        auto p = splitKLeft(root->left, k); root->left = p.second;
        recalc(root); return {p.first, root};
    }
}
Node* findByKey(Node* root, pkey_t key) {
    if (root==nullptr) return nullptr; push(root); if (root->pk==key) return root;
    return key < root->pk ? findByKey(root->left, key) : findByKey(root->right, key);
}
Node* findByPos(Node* root, int pos) {
    if (root==nullptr) return nullptr; int leftSz = size(root->left); push(root); if (pos==leftSz) return root;
    return pos < leftSz ? findByPos(root->left, pos) : findByPos(root->right, pos-leftSz-1);
}
//usable both in implicit and normal treap (for normal all left keys < all right keys)
Node* merge(Node* root1, Node* root2) {
    if (!root1) return root2; if (!root2) return root1;
    push(root1); push(root2);
    if (root1->sk>root2->sk) {
        root1->right = merge(root1->right, root2); recalc(root1); return root1;
    } else {
        root2->left = merge(root1, root2->left); recalc(root2); return root2;
    }//for persist, if sz1=L, sz2=R, choose root 1 with prob L/(L+R). Otherwise O(N) contertest
}
Node* getLeftMost(Node* root) {if (root->left) return getLeftMost(root->left); else return root;}
//Optimized typical ops (for NOT implicit treap)
Node* insert(Node* root, Node* node) {//TRIVIAL WAY: merge(<x, merge(newItem, >=x))
    if (!root) return node; push(root);
    if (node->sk < root->sk) {
        auto res = split(root, node->pk); node->l = res.first; node->r = res.second; return node;
    } else { return insert(node->pk < root->pk ? node->l : node->r, it);}
}
Node* erase(Node* root, pkey_t key) {//TRIVIAL WAY: split(x), merge(rootR->l, rootR->r), merge back
    assert(root); push(root);
    if (root->pk == key) return merge(root->left, root->right);
    else erase(key < root->pk ? root->left : root->right, key);
}
//to apply group op: cut tree with segment and put change to root. MAGIC works
//to make rq - same thing (
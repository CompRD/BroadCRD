//
// not well tested (or probably thought out) tree class
//
// this made some specficic experiments easy for me, but you might be better
// off just using a Digraph
//


template <class T>
class Tree {
public:
    Tree() = delete;
    Tree( T const& node ) : mNode(node), mLeft(nullptr), mRight(nullptr),
            mParent(nullptr) {};
    Tree* Left() const { return mLeft; }
    Tree* Right() const { return mRight; }
    T const& Node() const { return mNode; }
    void AddLeft( T const& node ) {
        delete mLeft;
        mLeft = new Tree(node);
        mLeft->mParent = this;
    }
    void AddRight( T const& node ) {
        delete mRight;
        mRight = new Tree(node);
        mRight->mParent = this;
    }
    enum class Order { FIRST, MIDDLE, LAST };
    String printOrder(Order o) {
        switch (o) {
        case Order::FIRST: return "FIRST"; break;
        case Order::MIDDLE: return "MIDDLE"; break;
        case Order::LAST: return "LAST"; break;
        }
        return "BUG";
    }
    Tree* DFWalk(std::function<bool(Tree*)> visit, Order o ) {
        std::vector<Tree*> stack;
        // mark where we've trodden.  we could have a set here instead, but i'm
        // thinking about storing something else... the value's never used right
        // now.
        std::map<Tree*, bool> explored;
        stack.push_back(this);
        while ( stack.size() ) {
            auto node = stack.back();
            stack.pop_back();
            ForceAssert(node);
            if ( explored.find( node ) != explored.end() ) {
                if ( visit(node ) ) return node;
            } else {
                explored[node]=true;
                if ( o == Order::LAST ) stack.push_back(node);
                if ( node->mRight && explored.find(node->mRight) == explored.end() )
                    stack.push_back(node->mRight);
                if ( o == Order::MIDDLE ) stack.push_back(node);
                if ( node->mLeft && explored.find(node->mRight) == explored.end() )
                    stack.push_back(node->mLeft);
                if ( o == Order::FIRST ) stack.push_back(node);
            }
        }
        return nullptr;
    }

    // really just for testing
    Tree* RecursiveDFWalk(std::function<bool(Tree*)> visit, Order o ) {
        if ( o == Order::FIRST) { if ( visit(this) ) return this; }
        if ( this->mLeft ) {
            auto tmp2 = this->mLeft->RecursiveDFWalk(visit, o);
            if ( tmp2 ) return tmp2;
        }
        if ( o == Order::MIDDLE ) { if ( visit(this) ) return this; }
        if ( this->mRight ) {
            auto tmp3 = this->mRight->RecursiveDFWalk(visit, o);
            if ( tmp3 ) return tmp3;
        }
        if ( o == Order::LAST ) { if ( visit(this) ) return this; }
        return nullptr;
    }

    void DeleteSubree() {
        DFWalk( [](Tree* node) { delete node; },  Order::LAST );
    }

    Tree* BFWalk( std::function<bool(Tree*)> visit ) {
        std::queue<Tree*> queue;
        queue.push(this);
        while ( queue.size() ) {
            auto node = queue.front();
            queue.pop();
            ForceAssert(node);
            if ( node->mLeft ) queue.push(node->mLeft);
            if ( node->mRight ) queue.push(node->mRight);
            if ( visit(node) ) return node;
        }
        return nullptr;
    }

    Tree* BinaryWalk( std::function<int(Tree*)> compare3way ) {
        // compare3way function returns -1 meaning "go left"
        // or +1 meaning "go right"
        // or 0 meaning "stop"
        Tree* current = this;
        while (current != nullptr) {
            int res = compare3way( current );
            if ( res == -1 ) current = current->Left();
            else if ( res == 1 ) current = current->Right();
            else if ( res == 0 ) break;
            else FatalErr("BinaryWalk comparator gave a bad result: "+ToString(res));
        }
        return current;
    }


    void SetLeft(Tree* const left) { mLeft = left; };
    void SetRight(Tree* const right) { mRight = right; };
    bool isLeft() const { return mLeft != nullptr; }
    bool isRight() const { return mRight != nullptr; }

private:
    T mNode;
    Tree* mLeft;
    Tree* mRight;
    Tree* mParent;
};

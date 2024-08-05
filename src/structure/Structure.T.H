template<class Key,class Value>
BinaryTree<Key,Value>& Foam::BinaryTree<Key,Value>::operator=
(
    const BinaryTree& rhs
)
{
    root = std::make_unique<Node>(*(rhs.root));
    return *this;
}

template<class Key,class Value>
Foam::BinaryTree<Key,Value>::Node::Node
(
    const Node& other
)
{
    key = other.key;
    value = other.value;
    if(other.leftChild)
        leftChild = std::make_unique<Node>(*(other.leftChild));
    if(other.rightChild)
        rightChild = std::make_unique<Node>(*(other.rightChild));
}

template<class Key,class Value>
const Value* Foam::BinaryTree<Key,Value>::find
(
    const Key& key
)
{
    auto recursiveGoDown = [&](Node* curr)
    {
        Node* child;
        if(key<=curr->key)
            child = curr->leftChild.get();
        else
            child = curr->rightChild.get();
        
        if(child!=nullptr)
            return recursiveGoDown(child);
        else
            return curr;
    };
    
    Value* value = nullptr;
    if(root)
    {
        Node* curr = recursiveGoDown(curr);
        value = curr->value;
    }
    return value;
}

template<typename T>
std::unique_ptr<List<List<T>>> Foam::Structure::broadcastHaloFields
(
    const GeometricField<T,fvPatchField,volMesh>& fieldData
)
{
    auto haloFieldPtr = std::make_unique<List<List<T>>>(Pstream::nProcs());
    List<List<T>>& haloField = *haloFieldPtr;
    List<T>& ownHaloField = haloField[Pstream::myProcNo()];
    
    const DynamicList<CellDescription>& ownHaloCells = getHaloCellList(Pstream::myProcNo());    
    ownHaloField.setSize(ownHaloCells.size());
    
    for(label haloCellInd=0; haloCellInd<ownHaloCells.size(); haloCellInd++)
    {
        ownHaloField[haloCellInd] = fieldData[ownHaloCells[haloCellInd].index];
    }
    
    return haloFieldPtr;
}

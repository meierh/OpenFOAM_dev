#include "MeshTree.H"

Foam::MeshTree::MeshTree
(
    const fvMesh& mesh
):
mesh(mesh)
{
    computeCellBBs();
}

void Foam::MeshTree::computeCellBBs()
{
    cellBBs.clear();
    const cellList& cells = mesh.cells();
    for(label cellInd=0; cellInd<mesh.cells().size(); cellInd++)
    {
        const cell& oneCell = cells[cellInd];
        BoundingBox bb(mesh,oneCell);
        cellBBs.append({cellInd,bb,{nullptr}});
    }
}

void Foam::MeshTree::constructTree()
{
    root = std::make_unique<node>();
    for(CellData& cData : cellBBs)
    {
        root->data.append(&cData);
    }    
    std::vector<node*> front = {root.get()};
    bool contSplit = true;
    while(contSplit)
    {
        std::vector<node*> newFront;
        for(node* n : front)
        {
            OctoSplit<bool> splitted = splitOutNode(n);
            for(label dx=0; dx<2; dx++)
            {
                for(label dy=0; dy<2; dy++)
                {         
                    for(label dz=0; dz<2; dz++)
                    {
                        if(splitted[dx][dy][dz])
                            newFront.push_back(n->leafs[dx][dy][dz].get());
                    }
                }
            }
        }
        if(newFront.size()==0)
            contSplit = false;
        front = newFront;
    }
}

Foam::BoundingBox Foam::MeshTree::minMaxSpan
(
    const DynamicList<CellData*>& boxes
) const
{
    if(boxes.size()==0)
        FatalErrorInFunction<<"Boxes size is zero"<<exit(FatalError);
    
    vector smaller = std::get<1>(*(boxes[0])).getSmaller();
    vector larger = std::get<1>(*(boxes[0])).getLarger();
    for(const CellData* cData : boxes)
    {
        for(label dim=0; dim<3; dim++)
        {
            smaller[dim] = std::min(smaller[dim],std::get<1>(*cData).getSmaller()[dim]);
            larger[dim] = std::min(larger[dim],std::get<1>(*cData).getLarger()[dim]);
        }
    }
    return BoundingBox(smaller,larger);
}

Foam::MeshTree::OctoSplit<bool> Foam::MeshTree::splitOutNode
(
    node* oneNode
)
{
    oneNode->box = minMaxSpan(oneNode->data);
    oneNode->split = 0.5*(oneNode->box.getSmaller()+oneNode->box.getLarger());
    
    for(label dx=0; dx<2; dx++)
    {
        for(label dy=0; dy<2; dy++)
        {
            for(label dz=0; dz<2; dz++)
            {
                oneNode->leafs[dx][dy][dz] = std::make_unique<node>();
            }
        }
    }    
    
    for(CellData* cData : oneNode->data)
    {
        labelList cDataSideX = side(cData,oneNode->split[0],0);
        labelList cDataSideY = side(cData,oneNode->split[1],1);
        labelList cDataSideZ = side(cData,oneNode->split[2],2);

        for(label dx : cDataSideX)
        {
            for(label dy : cDataSideY)
            {
                for(label dz : cDataSideZ)
                {
                    oneNode->leafs[dx][dy][dz]->data.append(cData);
                }
            }
        }
    }
    
    OctoSplit<bool> splitted;
    for(label dx=0; dx<2; dx++)
    {
        for(label dy=0; dy<2; dy++)
        {         
            for(label dz=0; dz<2; dz++)
            {
                if(oneNode->leafs[dx][dy][dz]->data.size() != oneNode->data.size())
                    splitted[dx][dy][dz] = true;
                else
                {
                    for(CellData* cData : oneNode->leafs[dx][dy][dz]->data)
                    {
                        std::get<2>(*cData).push_back(oneNode->leafs[dx][dy][dz].get());
                    }
                }
            }
        }
    }
    oneNode->data.clear();
    
    return splitted;
}

Foam::labelList Foam::MeshTree::side
(
    const CellData* cData,
    scalar value,
    label dim
) const
{
    const vector& smaller = std::get<1>(*cData).getSmaller();
    const vector& larger = std::get<1>(*cData).getLarger();
    if(value < smaller[dim])
        return {0};
    else if(larger[dim] < value)
        return {1};
    else
        return {0,1};
}

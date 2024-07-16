#include "LineStructure.H"

Foam::LineStructureStatus::LineStructureStatus()
{
    auto setFollows = [&](Node& method)
    {
        for(Node* methodPtr : method.dependencies)
            methodPtr->following.push_back(&method);
    };
    setFollows(initialPoints);
    setFollows(markers);
    setFollows(markersRefined);
    setFollows(markersVolume);
    setFollows(markersReduction);
    setFollows(markerMesh);
    setFollows(markersCollected);
    setFollows(markersCellWeight);
    setFollows(markersHaloCollect);
    setFollows(markersHaloExchange);
    setFollows(markersWeight);            
}

void Foam::LineStructureStatus::execValid
(
    const Node& method
)
{
    for(Node* depend : method.dependencies)
        if(!(depend->value))
        {
            Pout<<"Missing:"<<depend->name<<" from:"<<method.name<<Foam::endl;
            FatalErrorInFunction<<"Not all dependencies matched!"<<exit(FatalError);
        }
}

void Foam::LineStructureStatus::executed
(
    Node& method
)
{
    std::function<void(Node&)> recursiveUndo;
    recursiveUndo = [&](Node& method)
    {
        method.value = false;
        for(Node* methodPtr : method.following)
            recursiveUndo(*methodPtr);
    };
    recursiveUndo(method);
    method.value = true;
}

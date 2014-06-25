//
//  main.cpp
//  HW4
//

#include <iostream>
#include <fstream>
#include <ctime>
#include <random>
#include <vector>
#include <map>
#include <string.h>
#include <utility>
#include <functional>
#include <limits>
#include <queue>


//Note this code includes previous assignments (HW2 and HW3) The new additions are marked with HW4 comment line

enum Side {
    NORTH = -1,
    SOUTH = -2,
    WEST = -3,
    EAST = -4
};

//------------------------------------------------------------------------------
// EDGE STRUCTURE
//------------------------------------------------------------------------------
// Edge structure to store edge info
struct Edge {
    int originNode;
    int destinationNode;
    bool belongsToPlayerOne;
    bool belongsToPlayerTwo;
    double cost;
};

//------------------------------------------------------------------------------
// OPERATOR OVERLOADING << FOR TYPE EDGE
//------------------------------------------------------------------------------
// operator overloding of operator << to provide convinient output of type Edge
std::ostream& operator<<(std::ostream& os, const Edge& edge)
{
    // write obj to stream
    std::cout << "o: " << edge.originNode << " d: " << edge.destinationNode << " cost: " << edge.cost << "\n";
    return os;
}


//------------------------------------------------------------------------------
// GRAPH CLASS
//------------------------------------------------------------------------------
// class graph that implements the creation of random graph as well as shortest path algorithm
class HEXGraph
{
private:
    double density; // density of grpah
    int numberOfVertices; // number of nodes in the graph
    int numberOfEdges; // number of edges in the graph
    int numberOfSideNodes;
    std::map<int, std::vector <Edge>> adjacencyList; // map with key representing node id and vector representing node's connections to toher nodes with their distances
public:
    std::map<int, std::vector <Edge>>& getAdjacencyList() {return adjacencyList;} // returns adjacencyList of the graph
    int getNumberOfVertices() {return numberOfVertices;} // returns the number of vertices in the graph
    int getNumberOfEdges() {return numberOfEdges;} // returns the number of edges in the graph
    int getNumberOfSideNodes() {return numberOfSideNodes;}
    void addEdge(int originNode, int destinationNode, double edgeCost, bool belongsToPlayerOne = false, bool belongsToPlayerTwo = false); //adds to G the edge origin to destiantion with given edgeCost
    HEXGraph(int numberOfSideNodes = 7); //default constructor
    bool playerWin(int player);
    void markEdges(int horizontal, int vertical, int player);
};




//------------------------------------------------------------------------------
// OPERATOR OVERLOADING << FOR TYPE GRAPH
//------------------------------------------------------------------------------
// operator overloding of operator << to provide convinient output of type Graph
std::ostream& operator<<(std::ostream& os, HEXGraph& graph)
{
    // write obj to stream
    // using iterator object to loop through the map
    std::map<int, std::vector <Edge>>::iterator it1;
    for (it1 = graph.getAdjacencyList().begin(); it1 != graph.getAdjacencyList().end(); ++it1)
    {
        // using iterator object to loop through each vector
        for(std::vector<Edge>::iterator it2 = (it1->second).begin(); it2 != (it1->second).end(); ++it2)
        {
            std::cout << (*it2); //output type Edge
        }
        
        // Make sure you don't modify table here or the iterators will not work as you expect
    }
    std::cout << "number of edges: " << graph.getNumberOfEdges() << "\n";
    return os;
}

//------------------------------------------------------------------------------
// NODE STRUCTURE
//------------------------------------------------------------------------------
// Node structure for the Dijkstra algorithm
struct Node {
    int id;
    double distanceToSource; // Unknown distance from source to this node
    bool visited;        // Nodes have not been visited
    int previousNode;    // Previous node in optimal path from source
    int belongsToPlayer; //will show who made the move on this node
};

//------------------------------------------------------------------------------
// OPERATOR OVERLOADING << FOR TYPE EDGE
//------------------------------------------------------------------------------
// operator overloding of operator << to provide convinient output of type Edge
std::ostream& operator<<(std::ostream& os, const Node& node)
{
    // write obj to stream
    std::cout << node.previousNode << " \t-> " << node.id << " \tcost: " << node.distanceToSource << "\n";
    return os;
}

//------------------------------------------------------------------------------
// SUPPORTING STRUCT FOR DIJKSTRA ALGORITHM PRIORITY QUEUE
//------------------------------------------------------------------------------
// comparison operator used in the priority_queue in Graph::findShortestPath
struct compare
{
    bool operator()(const Node& l, const Node& r)
    {
        return l.distanceToSource > r.distanceToSource;
    }
};



// HW4
// HEXGraph should be a sublass of graph with additional properties and methods like
// - init with hex structure of size N
// - mark node i as Blue or Red
// - method to check if hex graph has a branch from one edge to another for particular color

/* HW 4 expectations:

1. Be able to draw the board using ASCII symbols and a given size, such as 7 by 7 or 11 by 11.
2. Input a move and determine if a move is legal.
3. Determine who won.

 X — . — . — . — .
  \ / \ / \ / \ / \
   . — . — . — . — .

*/

//------------------------------------------------------------------------------
// GRAPH::ADDEDGE
//------------------------------------------------------------------------------
// method to add undirectional edge to the graph's adjacency list
void HEXGraph::addEdge(int originNode, int destinationNode, double edgeCost, bool belongsToPlayerOne, bool belongsToPlayerTwo)
{
    // TODO: first check if the edge is not added already
    // init edge
    Edge e;
    e.originNode = originNode;
    e.destinationNode = destinationNode;
    e.cost = edgeCost;
    e.belongsToPlayerOne = belongsToPlayerOne;
    e.belongsToPlayerTwo = belongsToPlayerTwo;
    // add edge to adjacencyList
    this->adjacencyList[originNode].push_back(e);
    
    // enter the same edge but for opposite direction (because it is undirectional graph)
    e.originNode = destinationNode;
    e.destinationNode = originNode;
    //e.cost = edgeCost;
    this->adjacencyList[destinationNode].push_back(e);
    
    // increarse number of edges in the graph by 2 (because we just added 2 edges)
    this->numberOfEdges+=2;
}

//marks al edges of the node  as belonging to the player (to find out if he won or not)
void HEXGraph::markEdges(int horizontal, int vertical, int player)
{
    int node = horizontal + vertical * this->numberOfSideNodes;
    // for each neighbor node of current node:
    for(std::vector<Edge>::iterator it = this->adjacencyList[node].begin(); it != this->adjacencyList[node].end(); ++it)
    {
        if (player == 1) it->belongsToPlayerOne = true;
        else it->belongsToPlayerTwo = true;
    }
}

//------------------------------------------------------------------------------
// HEXGRAPH::HEXGRPAH CONSTRUCTOR (GENRATES EMPTY HEX GRAPH WITH number of nodes on one side)
//------------------------------------------------------------------------------
// constructor of class graph that creates random graph with paameters:
// int numberOfVertices - number of nodes in the grpah
// double density - density ofthe graph
// double rangeStart - start of range for the edge cost
// double rangeEnd - end of range for the edge cost
HEXGraph::HEXGraph(int numberOfSideNodes)
{
    // init
    this->numberOfEdges = 0;
    this->numberOfSideNodes = numberOfSideNodes;
    // set number of vertices
    this->numberOfVertices = numberOfSideNodes * numberOfSideNodes; //calculate number of nodes
    
    // looping through each side node per one row
    for (int i = 0; i < numberOfVertices; i++)
    {
        // looping through each side node in that row for ech column
        for (int j = i+1; j < numberOfVertices; j++)
        {
            // general connections: (i, j-1),(i,j+1),(i+1,j-1),(i+1,j),(i-1,j-1),(i-1,j)
            //to convert i or j to board coords need to modulo of dividing by numberOfSideNodes
            int n = i % numberOfSideNodes;
            int m = j % numberOfSideNodes;
            
            //check conditions for neighbor nodes
            if  (((i == j-1) && (m != 0)) || (i == j-numberOfSideNodes) || ((i == j-numberOfSideNodes+1) && (n != 0)))
                this->addEdge(i, j, 1);
        }
        
        //virtual nodes
        int n = i % numberOfSideNodes;
        
        if (n == 0) this->addEdge(WEST, i, 1, true, false); // virtual node for west side
        if (n == numberOfSideNodes-1) this->addEdge(EAST, i, 1, true, false); // virtual node for east side
        if (i < numberOfSideNodes) this->addEdge(NORTH, i, 1, false, true); // virtual node for north side
        if (i > numberOfVertices - numberOfSideNodes) this->addEdge(SOUTH, i, 1, false, true); // virtual node for south side
    }
}

//------------------------------------------------------------------------------
// GRAPH::FIND SHORTEST PATH (DIJKSTRA ALGORITHM)
//------------------------------------------------------------------------------
// class implementing Dijkstra algorithm
// for this node find shortest path to another node if exists
bool HEXGraph::playerWin(int player)
{
    std::vector<Node> nodes; // this vector will store all nodes of the graph using structure Node (see above)
    // initializations
    for (int i = 0; i < this->getNumberOfVertices(); i++)
    {
        Node n;
        n.id = i; // current node id equal i
        n.distanceToSource = std::numeric_limits<double>::infinity();
        n.visited = false;
        n.previousNode = -1.0;
        n.belongsToPlayer = 0;
        nodes.push_back(n); // add node n to the nodes vector
    }
    int source, destination;
    if (player == 1)
    {
        source = EAST;
        destination = WEST;
    }
    else if (player == 2)
    {
        source = NORTH;
        destination = SOUTH;
    }
    else
    {
        return false;
        std::cout << "no such player";
    }
    nodes[source].distanceToSource = 0;
    
    std::priority_queue<Node, std::vector<Node>, compare> priorityQueue; // priority queue that always have the node with the shortet distnce to the source on the top (thanks to the compare operator (see above this method)
    priorityQueue.push(nodes[source]); // add current node which is the source node
    
    // main loop
    while (!priorityQueue.empty())
    {
        Node currentNode;
        currentNode = priorityQueue.top(); // curnet node := vertex in priority queue with smallest distance in dist[] and has not been visited;  // Source node in first case
        
        if (currentNode.id == destination) // if the curent node is the destination node the shortest path is found
        {
            return true; // nodes[destination].distanceToSource; // return shortest path
        }
        
        priorityQueue.pop(); // remove currentNode from priority queue;
        nodes[currentNode.id].visited = true; // mark this node as visited
        
        // for each neighbor node of current node:
        for(std::vector<Edge>::iterator it = this->adjacencyList[currentNode.id].begin(); it != this->adjacencyList[currentNode.id].end(); ++it)
        {
            double currentDistance = nodes[currentNode.id].distanceToSource + (*it).cost; // accumulate shortest dist from source
            // if current disrance is below distance to source for current destiantoin node and current destination node was not visited then
            if ((currentDistance < nodes[(*it).destinationNode].distanceToSource) && (!nodes[(*it).destinationNode].visited) && (nodes[(*it).destinationNode].belongsToPlayer = player))
            {
                nodes[(*it).destinationNode].distanceToSource = currentDistance; // set the shortest dist from src to current destiantion node
                nodes[(*it).destinationNode].previousNode = currentNode.id; // set previuos node of current destination node to current node =)
                priorityQueue.push(nodes[(*it).destinationNode]); // Add unvisited curent destiantion node into the priority queue to be process
            } //end if
        }   //end for
        
    } //end while
    
    return false; //nodes[destination].distanceToSource; // return distance to source
}



//------------------------------------------------------------------------------
// MAIN
//------------------------------------------------------------------------------
int main(int argc, const char * argv[])
{

    int i, horizontal, vertical;
    
    while (true)
    {
        std::cout << "Please enter number of side nodes\n";
        std::cin >> i;
        HEXGraph hexGraph(i);
        std::cout << hexGraph;
        std::cout << "You are player number 1 (EAST -> WEST)\n";
        std::cout << "Please enter your move as 2 integers: horizontal and vertical\n";
        std::cin >> horizontal >> vertical;
        hexGraph.markEdges(horizontal, vertical, 1);
        if (hexGraph.playerWin(1)) std::cout << "You won!!!";
    }
    
    return 0;
}


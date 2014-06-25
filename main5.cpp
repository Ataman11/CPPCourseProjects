//
//  main.cpp
//  HW5
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
#include <unordered_map>


enum Side {
    WEST = 0,
    EAST = 1,
    NORTH = 2,
    SOUTH = 3
};

//------------------------------------------------------------------------------
// EDGE STRUCTURE
//------------------------------------------------------------------------------
// Edge structure to store edge info
struct Edge {
    int originNode;
    int destinationNode;
    int belongsToPlayer;
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
// OPERATOR OVERLOADING << FOR TYPE NODE
//------------------------------------------------------------------------------
// operator overloding of operator << to provide convinient output of type Edge
std::ostream& operator<<(std::ostream& os, const Node& node)
{
    // write obj to stream
    std::cout << node.previousNode << " \t-> " << node.id << " \tcost: " << node.distanceToSource << "\n";
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
    std::unordered_map<int, std::vector <Edge>> adjacencyList; // map with key representing node id and vector representing node's connections to toher nodes with their distances
    std::vector<Node> nodesList; //contains nodes to determine if the nodes belongs to player
    std::vector<int> shuffledNodes; //helper property for AI
public:
    std::unordered_map<int, std::vector <Edge>>& getAdjacencyList() {return adjacencyList;} // returns adjacencyList of the graph
    int getNumberOfVertices() {return numberOfVertices;} // returns the number of vertices in the graph
    int getNumberOfEdges() {return numberOfEdges;} // returns the number of edges in the graph
    int getNumberOfSideNodes() {return numberOfSideNodes;}
    void addEdge(int originNode, int destinationNode, double edgeCost, int belongsToPlayer = 0); //adds to G the edge origin to destiantion with given edgeCost
    void addNode(int node, int belongsToPlayer = 0);
    HEXGraph(int numberOfSideNodes = 7); //default constructor
    HEXGraph(HEXGraph const &graph); //copy constructor
    bool playerWin(int player);
    bool markEdges(std::string moveString, int player);
    void markNode(int node, int player);
    void print();
    void computerMove();
};




//------------------------------------------------------------------------------
// OPERATOR OVERLOADING << FOR TYPE GRAPH
//------------------------------------------------------------------------------
// operator overloding of operator << to provide convinient output of type Graph
std::ostream& operator<<(std::ostream& os, HEXGraph& graph)
{
    // write obj to stream
    // using iterator object to loop through the map
    std::unordered_map<int, std::vector <Edge>>::iterator it1;
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

//helper function
std::string playerToStr(int palyer)
{
    std::string str;
    switch (palyer) {
        case 1:
            str = "X";
            break;
        case 2:
            str = "O";
            break;
        default:
            str = ".";
            break;
    }
    return str;
}

//hex board print method
void HEXGraph::print()
{
    char rowName = 'A';
    
    std::cout << "\n";
    
    std::cout << "  ";
    for (int i = 0; i < this->numberOfSideNodes; i++)
    {
        std::cout << i << "   ";
    }
    std::cout << "\n";
    
    int node = 0;
    
    // looping through each side node per one row to draw the graph
    for (int i = 0; i < this->numberOfSideNodes-1; i++)
    {
        for (int j = 0; j < i*2; j++)
        {
            std::cout << " ";
        }
        
        std::cout << rowName << " " << playerToStr(this->nodesList[node].belongsToPlayer);
        node++;
        // looping through each side node in that row for each column
        for (int j = 0; j < numberOfSideNodes-1; j++)
        {
            std::cout << " — " << playerToStr(this->nodesList[node].belongsToPlayer);
            node++;
        }
        std::cout << "\n";
        
        for (int j = 0; j < i*2+2; j++)
        {
            std::cout << " ";
        }
        std::cout << " \\";
        // looping through each side node in that row for each column
        for (int j = 0; j < numberOfSideNodes-1; j++)
        {
            std::cout << " / \\";
        }
        std::cout << "\n";
        
        rowName++;
    }
    
    for (int j = 0; j < (numberOfSideNodes - 1)*2; j++)
    {
        std::cout << " ";
    }
    std::cout << rowName << " " << playerToStr(this->nodesList[node].belongsToPlayer);
    node++;
    // looping through each side node in that row for each column
    for (int j = 0; j < numberOfSideNodes-1; j++)
    {
        std::cout << " — " << playerToStr(this->nodesList[node].belongsToPlayer);
        node++;
    }
    std::cout << "\n";
    std::cout << "\n";
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

//------------------------------------------------------------------------------
// GRAPH::ADDEDGE
//------------------------------------------------------------------------------
// method to add undirectional edge to the graph's adjacency list
void HEXGraph::addEdge(int originNode, int destinationNode, double edgeCost, int belongsToPlayer)
{
    // TODO: first check if the edge is not added already
    // init edge
    Edge e;
    e.originNode = originNode;
    e.destinationNode = destinationNode;
    e.cost = edgeCost;
    e.belongsToPlayer = belongsToPlayer;
    e.belongsToPlayerOne = false;
    e.belongsToPlayerTwo = false;
    if (belongsToPlayer == 1) e.belongsToPlayerOne = true;
    else if (belongsToPlayer == 2) e.belongsToPlayerTwo = true;
    
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

//------------------------------------------------------------------------------
// GRAPH::ADDNODE
//------------------------------------------------------------------------------
// method to add undirectional edge to the graph's adjacency list
void HEXGraph::addNode(int node, int belongsToPlayer)
{
    // TODO: first check if the edge is not added already
    // init edge
    Node n;
    n.id = node;
    n.belongsToPlayer = belongsToPlayer;

    this->nodesList.push_back(n);
}




//marks al edges of the node  as belonging to the player (to find out if he won or not)
bool HEXGraph::markEdges(const std::string moveString, int player)
{
    if (moveString.length() < 2)
    {
        std::cout << "!!! Please enter move in proper format.\n";
        return false;
    }
    char c1 = moveString[0];
    char c2 = moveString[1];
    
    char maxChar = static_cast<char>(this->numberOfSideNodes + 65 - 1);
    //char maxNumber = static_cast<char>(this->numberOfSideNodes + 49);
    //check if c1 is letter
    //if (c1 <= 'A') || (c1 >= maxChar) || (c1 <= '1') || (c1 >= maxNumber) ||
    
    int column, row;
    column = atoi(&c1);
    //if string starts with int then column number was entered first
    if (column > 0)
    {
        column = std::stoi(moveString);
        if (column >= this->numberOfSideNodes)
        {
            std::cout << "!!! Please enter correct column number.\n";
            return false;
        }
        if (column > 9) row = static_cast<int>(moveString[2]) - 65;
        else row = static_cast<int>(moveString[1]) - 65;
        if ((row < 0) || (row >= this->numberOfSideNodes))
        {
            std::cout << "!!! Please enter correct row letter.\n";
            return false;
        }
    }
    //the letter was entered first
    else
    {
        row = static_cast<int>(moveString[0]) - 65;
        if ((row < 0) || (row >= this->numberOfSideNodes))
        {
            std::cout << "!!! Please enter correct row letter.\n";
            return false;
        }
        column = std::stoi(moveString.substr(1, 2));
        if ((column >= this->numberOfSideNodes) || (column < 0))
        {
            std::cout << "!!! Please enter correct column number.\n";
            return false;
        }
        
    }
    
    int node = column + row * this->numberOfSideNodes; //get the absolute number of node
    
    //mark node first
    if (this->nodesList[node].belongsToPlayer != 0)
    {
        std::cout << "!!! This spot is occupied already!\n";
        return false;
    }
    
    this->nodesList[node].belongsToPlayer = 1;
    // for each neighbor node of current node edge is also considered occupied by the player
    for(std::vector<Edge>::iterator it = this->adjacencyList[node].begin(); it != this->adjacencyList[node].end(); ++it)
    {
        it->belongsToPlayerOne = true;
    }
    
    return true;

}

//------------------------------------------------------------------------------
// HEXGRAPH::HEXGRPAH CONSTRUCTOR (GENRATES EMPTY HEX GRAPH WITH number of nodes on one side)
//------------------------------------------------------------------------------

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
        
        if (n == 0) this->addEdge(numberOfVertices + WEST, i, 1, 1); // virtual node for west side
        if (n == numberOfSideNodes-1) this->addEdge(numberOfVertices + EAST, i, 1, 1); // virtual node for east side
        if (i < numberOfSideNodes) this->addEdge(numberOfVertices + NORTH, i, 1, 2); // virtual node for north side
        if (i >= numberOfVertices - numberOfSideNodes) this->addEdge(numberOfVertices + SOUTH, i, 1, 2); // virtual node for south side
    }
    
    //create nodes list
    for (int i = 0; i < numberOfVertices; i++) this->addNode(i);
    //add the virtual nodes to nodes list
    /*WEST = 0,
    EAST = 1,
    NORTH = 2,
    SOUTH = 3*/
    this->addNode(numberOfVertices + WEST, 1);
    this->addNode(numberOfVertices + EAST, 1);
    this->addNode(numberOfVertices + NORTH, 2);
    this->addNode(numberOfVertices + SOUTH, 2);
    
    //init shuffle nodes vector
    for (int i = 0; i < this->numberOfVertices; i++) {
        this->shuffledNodes.push_back(i);
    }

}

//copy constructor
HEXGraph::HEXGraph(HEXGraph const &graph)
{
    this->numberOfVertices = graph.numberOfVertices; // number of nodes in the graph
    this->numberOfEdges = graph.numberOfEdges; // number of edges in the graph
    this->numberOfSideNodes = graph.numberOfSideNodes;
    std::unordered_map<int, std::vector <Edge>> newAdjacencyList(graph.adjacencyList);
    this->adjacencyList = newAdjacencyList; // map with key representing node id and vector representing node's connections to toher nodes with their distances
    std::vector<Node> newNodesList(graph.nodesList);
    this->nodesList = newNodesList; //contains nodes to determine if the nodes belongs to player
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
    for (int i = 0; i < this->getNumberOfVertices() + 4; i++)
    {
        Node n;
        n.id = i; // current node id equal i
        n.distanceToSource = std::numeric_limits<double>::infinity();
        n.visited = false;
        n.previousNode = -1.0;
        n.belongsToPlayer = 0;
        nodes.push_back(n); // add node n to the nodes map
    }
    int source, destination;
    if (player == 1)
    {
        source = numberOfVertices + EAST;
        destination = numberOfVertices + WEST;
    }
    else if (player == 2)
    {
        source = numberOfVertices + NORTH;
        destination = numberOfVertices + SOUTH;
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
            if ((currentDistance < nodes[(*it).destinationNode].distanceToSource) && (!nodes[(*it).destinationNode].visited) &&
                (this->nodesList[(*it).destinationNode].belongsToPlayer == player) && (((*it).belongsToPlayerOne && (player == 1)) || ((*it).belongsToPlayerTwo && (player == 2))))
            {
                nodes[(*it).destinationNode].distanceToSource = currentDistance; // set the shortest dist from src to current destiantion node
                nodes[(*it).destinationNode].previousNode = currentNode.id; // set previuos node of current destination node to current node =)
                priorityQueue.push(nodes[(*it).destinationNode]); // Add unvisited curent destiantion node into the priority queue to be process
            } //end if
        }   //end for
        
    } //end while
    
    return false; //nodes[destination].distanceToSource; // return distance to source
}

void HEXGraph::markNode(int node, int player)
{
    this->nodesList[node].belongsToPlayer = player;
    for(std::vector<Edge>::iterator it = this->adjacencyList[node].begin(); it != this->adjacencyList[node].end(); ++it)
    {
        if (player == 1)  it->belongsToPlayerOne = true;
        else it->belongsToPlayerTwo = true;
    }
}


//AI - computer move algorithm

void HEXGraph::computerMove()
{
    //dummy random move generator is commented out
    //left it just in case for debugging puposes
    /*while (true)
    {
        int row = rand() % this->numberOfSideNodes;
        int column = rand() % this->numberOfSideNodes;
        
        int node = column + row * this->numberOfSideNodes; //get the absolute number of node
        
        //mark node first
        if (this->nodesList[node].belongsToPlayer == 0)
        {
            this->nodesList[node].belongsToPlayer = 2;
            // for each neighbor node of current node:
            for(std::vector<Edge>::iterator it = this->adjacencyList[node].begin(); it != this->adjacencyList[node].end(); ++it)
            {
                it->belongsToPlayerTwo = true;
            }
            break;
        }
    }*/
    
    int numberOfWinsMax = 0;
    int nunberOfWins;
    int nodeMax;
    //for (all unoccupied spaces on the board i)
    for(std::vector<Node>::iterator it = this->nodesList.begin(); it != this->nodesList.end(); ++it)
    //for (int i = 0; i++; i < this->nodesList.size())
    {
        if ((*it).belongsToPlayer == 0) //we check if node is unnocupied
        {
            nunberOfWins = 0;
            //for (all Monte Carlo trials)
            for (int j = 0; j < 1000; ++j)
            {
                //-make a copy of your board graph (at the state it is in before the computer makes its move)
                HEXGraph copyOfGraph(*this); //vector<int> new_(original); - fast copy
                //-in this copy
                //- update the board giving the computer's player unoccupied space i
                copyOfGraph.nodesList[(*it).id].belongsToPlayer = 2;
                for(std::vector<Edge>::iterator it2 = copyOfGraph.adjacencyList[(*it).id].begin(); it2 != copyOfGraph.adjacencyList[(*it).id].end(); ++it2)
                {
                    it2->belongsToPlayerTwo = true;
                }
                //- randomly fill the remaining unoccupied spaces (giving the human the next move)
                
                /*for(std::vector<Node>::iterator it3 = copyOfGraph.nodesList.begin(); it3 != copyOfGraph.nodesList.end(); ++it3)
                {
                    if ((*it3).belongsToPlayer == 0)
                    {
                        int randomPlayer = rand() % 2 + 1;
                        copyOfGraph.nodesList[(*it3).id].belongsToPlayer = randomPlayer;
                        for(std::vector<Edge>::iterator it4 = copyOfGraph.adjacencyList[(*it3).id].begin(); it4 != copyOfGraph.adjacencyList[(*it3).id].end(); ++it4)
                        {
                            if (randomPlayer == 1) it4->belongsToPlayerOne = true;
                            else it4->belongsToPlayerTwo = true;
                        }
                    }
                }*/
                
                //improvement - only randomize half of the moves and check if palyer 2 won.
                std::random_shuffle(this->shuffledNodes.begin(),this->shuffledNodes.end());
                
                for (int k = 0; k < this->numberOfVertices / 2; k++) {
                    if (nodesList[shuffledNodes[k]].belongsToPlayer == 0)
                    {
                        copyOfGraph.nodesList[shuffledNodes[k]].belongsToPlayer = 2;
                        for(std::vector<Edge>::iterator it4 = copyOfGraph.adjacencyList[shuffledNodes[k]].begin(); it4 != copyOfGraph.adjacencyList[shuffledNodes[k]].end(); ++it4)
                        {
                            it4->belongsToPlayerTwo = true;
                        }
                    }
                }
                
                
                //- determine who wins with as in HW4 (there is always a winner)
                if (copyOfGraph.playerWin(2))
                {
                    //- if the computer wins, nwins = nwins+1
                    nunberOfWins++;
                }
            }
            if (nunberOfWins > numberOfWinsMax)
            {
                numberOfWinsMax = nunberOfWins;
                nodeMax = (*it).id;
            }
        }
    }
    //Finally, the computer chooses the unoccupied space i = imax since it has the most MC wins of all the unoccupied spaces.
    this->nodesList[nodeMax].belongsToPlayer = 2;
    // for each neighbor node of current node:
    for(std::vector<Edge>::iterator it = this->adjacencyList[nodeMax].begin(); it != this->adjacencyList[nodeMax].end(); ++it)
    {
        it->belongsToPlayerTwo = true;
    }
}

//------------------------------------------------------------------------------
// MAIN
//------------------------------------------------------------------------------
int main(int argc, const char * argv[])
{

    /* initialize random seed: */
    srand (time(NULL));
    
    int i;
    std::string str;

    std::cout << "Please enter enter size of the board (example: 7)\n";
    //TODO: check the int
    std::cin >> i;
    HEXGraph hexGraph(i);
    //std::cout << hexGraph;
    std::cout << hexGraph;
    hexGraph.print();
    //char c[] = "33333";
    //std::cout << static_cast<int>(c);
    std::string str1 = "13dfsdf";
    std::cout << std::stoi(str1);
    
    std::cout << "You are player number 1 (EAST -> WEST)\n";
    std::cout << "Enter your move in format C5 or 1B.\n";
    
    while (true)
    {
        std::cout << "Please enter your move:\n";
        std::cin >> str;
        if (!hexGraph.markEdges(str, 1)) continue;
        
        if (hexGraph.playerWin(1))
        {
            std::cout << "************ You won!!! ***********\n";
            hexGraph.print();
            std::cout << "************ Game over ************\n";
            break;
        }
        hexGraph.computerMove();
        if (hexGraph.playerWin(2))
        {
            std::cout << "************ Computer won *********\n";
            hexGraph.print();
            std::cout << "************ Game over ************\n";
            break;
        }
        std::cout << "Computer moved:";
        hexGraph.print();
        
    }
    
    return 0;
}


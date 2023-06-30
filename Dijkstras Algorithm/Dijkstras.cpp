#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iomanip>

// Structure to represent an edge in the graph
struct Edge {
    int source;
    int destination;
    double distance;
};

// Class to represent a graph
class Graph {
private:
    int numVertices;
    std::vector<std::vector<Edge>> adjacencyList;

public:
    // Constructor
    Graph(int numVertices) : numVertices(numVertices), adjacencyList(numVertices) {}

    // Function to add an edge to the graph
    void addEdge(int source, int destination, double distance) {
        Edge edge;
        edge.source = source;
        edge.destination = destination;
        edge.distance = distance;
        adjacencyList[source].push_back(edge);
        std::swap(edge.source, edge.destination);
        adjacencyList[destination].push_back(edge);
    }

    // Function to generate a random graph with the given density and distance range
    void generateRandomGraph(double density, double distanceRange) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        std::uniform_real_distribution<> disDist(1.0, distanceRange);

        for (int i = 0; i < numVertices; ++i) {
            for (int j = i + 1; j < numVertices; ++j) {
                if (dis(gen) < density) {
                    double distance = disDist(gen);
                    addEdge(i, j, distance);
                }
            }
        }
    }

    // Function to calculate the shortest path between two vertices using Dijkstra's algorithm
    double shortestPath(int source, int destination) {
        std::vector<double> dist(numVertices, std::numeric_limits<double>::max());
        dist[source] = 0.0;

        std::vector<bool> visited(numVertices, false);

        for (int i = 0; i < numVertices - 1; ++i) {
            int minIndex = -1;
            double minDist = std::numeric_limits<double>::max();

            for (int v = 0; v < numVertices; ++v) {
                if (!visited[v] && dist[v] <= minDist) {
                    minIndex = v;
                    minDist = dist[v];
                }
            }

            visited[minIndex] = true;

            for (const auto& edge : adjacencyList[minIndex]) {
                double newDist = dist[minIndex] + edge.distance;
                if (!visited[edge.destination] && newDist < dist[edge.destination])
                    dist[edge.destination] = newDist;
            }
        }

        return dist[destination];
    }
};

// Function to calculate the average path length for a set of randomly generated graphs
double calculateAveragePathLength(int numVertices, double density, double distanceRange) {
    double totalPathLength = 0.0;
    int numPaths = 0;

    for (int destination = 1; destination <= numVertices; ++destination) {
        double pathLength = 0.0;
        int numValidPaths = 0;

        for (int source = 0; source < numVertices; ++source) {
            if (source == destination)
                continue;

            Graph graph(numVertices);
            graph.generateRandomGraph(density, distanceRange);

            double shortestPath = graph.shortestPath(source, destination);
            if (shortestPath != std::numeric_limits<double>::max()) {
                pathLength += shortestPath;
                numValidPaths++;
            }
        }

        if (numValidPaths > 0) {
            totalPathLength += (pathLength / numValidPaths);
            numPaths++;
        }
    }

    return totalPathLength / numPaths;
}

int main() {
    int numVertices = 50;
    double density1 = 0.2;
    double density2 = 0.4;
    double distanceRange = 10.0;

    double averagePathLength1 = calculateAveragePathLength(numVertices, density1, distanceRange);
    double averagePathLength2 = calculateAveragePathLength(numVertices, density2, distanceRange);

    std::cout << "Average Path Length (Density: 20%): " << std::fixed << std::setprecision(2) << averagePathLength1 << std::endl;
    std::cout << "Average Path Length (Density: 40%): " << std::fixed << std::setprecision(2) << averagePathLength2 << std::endl;

    return 0;
}

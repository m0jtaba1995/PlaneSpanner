/*
 * MIT License
 *
 * Copyright (c) 2024 Mojtaba Amani
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <random>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Point_2, EdgeWeightProperty> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::adjacency_iterator adjacency_iterator;

const double INF = std::numeric_limits<double>::infinity();

// Function to compute Euclidean distance between two points.
double euclidean_distance(const Point_2 &p1, const Point_2 &p2) {
    return CGAL::sqrt(CGAL::squared_distance(p1, p2));
}

// Function to compute Squared distance between two points for Rotating Calipers method.
double squared_distance(const Point_2 &p1, const Point_2 &p2) {
    return CGAL::squared_distance(p1, p2);
}

// Function to write generated plane t_spanner to a .txt file.
void write_spanner_to_file(const std::vector<std::pair<Point_2, Point_2>> &edges, std::string &file_name) {
    std::ofstream outfile(file_name);
    if (!outfile) {
        std::cerr << "Error: Unable to open file for writing: " << file_name << std::endl;
        return;
    }
    for (const auto &edge: edges) {
        outfile << edge.first.x() << " " << edge.first.y() << " " << edge.second.x() << " " << edge.second.y()
                << std::endl;
    }
    outfile.close();
}

// Function to find two neighbors of a point on a convex hull
std::pair<Point_2, Point_2> get_neighbors_on_convex_hull(const std::vector<Point_2> &points, const Point_2 &p) {
    std::pair<Point_2, Point_2> neighbors;

    // Find the index of point p in the convex hull
    auto it = std::find(points.begin(), points.end(), p);
    if (it != points.end()) {
        std::size_t index = std::distance(points.begin(), it);
        std::size_t size = points.size();

        // Find the clockwise neighbor
        neighbors.first = points[(index + 1) % size];

        // Find the counterclockwise neighbor
        neighbors.second = points[index == 0 ? size - 1 : index - 1];
    }

    return neighbors;
}

bool get_random_boolean_values() {
    // Seed the random number generator with a random device
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a uniform distribution for integers (0 or 1)
    std::uniform_int_distribution<> dist(0, 1);

    // Generate a random boolean value
    return dist(gen);
}

/*
 * Function to find the diameter of a polygon using the rotating calipers method in O(n log n) time
 * (in O(n) if points in convex position).
 * If your set of input points are not in convex position, set the isInConvexPosition = false so that the convex hull of input is find first.
 * Otherwise, set it true.
 */
std::pair<Point_2, Point_2> find_diametral_pair(const std::vector<Point_2> &points, bool is_points_in_convex_position) {
    std::vector<Point_2> convex_hull;
    if (!is_points_in_convex_position) {
        // Compute convex hull only if not in convex position
        CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convex_hull));
    } else {
        convex_hull = points; // Use points as convex hull
    }

    std::pair<Point_2, Point_2> max_distance_pair;
    double max_squared_distance = 0.0;
    std::size_t n = convex_hull.size();
    std::size_t j = 1;

    // Rotating calipers approach.
    for (std::size_t i = 0; i < n; i++) {
        while (true) {
            std::size_t next_j = (j + 1) % n;
            double current_squared_distance = squared_distance(convex_hull[i], convex_hull[j]);
            double next_squared_distance = squared_distance(convex_hull[i], convex_hull[next_j]);
            if (next_squared_distance > current_squared_distance) {
                j = next_j;
            } else {
                break;
            }
        }
        double current_squared_distance = squared_distance(convex_hull[i], convex_hull[j]);
        if (current_squared_distance > max_squared_distance) {
            max_squared_distance = current_squared_distance;
            max_distance_pair = std::make_pair(convex_hull[i], convex_hull[j]);
        }
    }
    return max_distance_pair;
}

/*
 * Function to find the approximate diametral pair of a set of points in convex position
 * (Janardan https://doi.org/10.1142/S021819599300021X)
 */
std::pair<Point_2, Point_2> find_approximate_diametral_pair(const std::vector<Point_2> &points, int c) {

    // maxPair will hold the pair of points with the maximum distance found.
    std::pair<Point_2, Point_2> maxPair;
    // maxDistance keeps track of the maximum distance found so far.
    double maxDistance = 0.0;

    // Iterate through c different coordinate systems.
    for (int i = 0; i < c; ++i) {
        // Calculate the rotation angle for the current coordinate system.
        double angle = i * M_PI / c;

        // Initialize minProj and maxProj with extreme values to find the minimum and maximum projections.
        double minProj = std::numeric_limits<double>::max();
        double maxProj = -std::numeric_limits<double>::max();

        // Initialize minPoint and maxPoint to store the points with min and max projections.
        Point_2 minPoint, maxPoint;

        // Iterate over all points to find the min and max projections in the current coordinate system.
        for (const auto &p: points) {
            // Compute the projection of point p along the axis defined by the current angle.
            double proj = p.x() * std::cos(angle) + p.y() * std::sin(angle);

            // Update minProj and minPoint if the current projection is smaller than the known minimum.
            if (proj < minProj) {
                minProj = proj;
                minPoint = p;
            }

            // Update maxProj and maxPoint if the current projection is larger than the known maximum.
            if (proj > maxProj) {
                maxProj = proj;
                maxPoint = p;
            }
        }

        // Calculate the Euclidean distance between minPoint and maxPoint.
        double distance = euclidean_distance(minPoint, maxPoint);

        // Update maxPair and maxDistance if the current distance is greater than the known maximum distance.
        if (distance > maxDistance) {
            maxDistance = distance;
            maxPair = {minPoint, maxPoint};
        }
    }
    return maxPair;
}

// Function to calculate the length of the convex hull.
double compute_perimeter(const std::vector<Point_2> &points) {
    double perimeter = 0.0;
    for (int i = 0; i < points.size(); i++) {
        perimeter += euclidean_distance(points[i], points[(i + 1) % points.size()]);
    }
    return perimeter;
}

// Function to find the proposed t_good point.
std::pair<double, Point_2> find_proposed_t_good(const std::vector<Point_2> &points, const double &perimeter) {
    std::size_t n = points.size();
    Point_2 goodPoint;
    double dilation = std::numeric_limits<double>::max();

    /*
     * A matrix to store the length of the path between two points A and B on the convex hull.
     * The main diagonal is zero.
     * All the entries above the main diagonal are the length of the path between A and B on the convex hull in the counterclockwise direction,
     * and the entries below the main diagonal are in the clockwise direction.
     */
    std::vector<std::vector<double>> distanceMatrix(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        distanceMatrix[i][i] = 0;
        for (int j = i + 1; j < n; ++j) {
            distanceMatrix[i][j] = distanceMatrix[i][j - 1] + euclidean_distance(points[j - 1], points[j]);
            distanceMatrix[j][i] = perimeter - distanceMatrix[i][j];
        }
    }
    for (std::size_t i = 0; i < n; i++) {
        double currentDilation = std::numeric_limits<double>::lowest();
        Point_2 currentPoint = points[i];

        for (std::size_t j = 0; j < n; j++) {
            if (j == i)
                continue;

            double t = std::min(distanceMatrix[i][j], distanceMatrix[j][i]) / euclidean_distance(points[i], points[j]);

            if (t >= currentDilation) {
                currentDilation = t;
                if (t > dilation)
                    break;
            }
        }

        if (currentDilation < dilation) {
            dilation = currentDilation;
            goodPoint = currentPoint;
        }
    }
    return std::make_pair(dilation, goodPoint);
}

// Function to compute the stretch factor for the input graph.
double compute_stretch_factor(const Graph &graph) {
    double maxStretch = 0.0;
    std::vector<std::vector<double>> distances(boost::num_vertices(graph),
                                               std::vector<double>(boost::num_vertices(graph), INF));
    /*
     * According to Boost document: (https://www.boost.org/doc/libs/)
     *  To find the shortest distance between every pair of vertices in the graph you have two choices:
     *  1- Use Floyd-Warshall algorithm for dense graphs. The time complexity is O(V^3).
     *  2- Use Johnson algorithm for sparse graphs. The time complexity is O(V E log V).
     *  These methods returns false if there is a negative weight cycle in the graph and true otherwise.
     *
     *  Both methods are available here and you can choose either one according to the type of graph you have.
     */
    //    boost::floyd_warshall_all_pairs_shortest_paths(graph, distances, boost::weight_map(get(boost::edge_weight, graph)));
    boost::johnson_all_pairs_shortest_paths(graph, distances, boost::weight_map(get(boost::edge_weight, graph)));

    /**
     * This part print the shortest path table for all pairs of vertices.
     *  Uncomment it if needed.
     */
    //    std::cout << "Shortest paths distance matrix:" << std::endl;
    //    for (auto &distance: distances) {
    //        for (double j: distance) {
    //            if (j == INF) {
    //                std::cout << std::setw(15) << "inf";
    //            } else {
    //                std::cout << std::setw(15) << j;
    //            }
    //        }
    //        std::cout << std::endl;
    //    }

    for (Vertex source = 0; source < boost::num_vertices(graph); source++) {
        for (Vertex target = 0; target < boost::num_vertices(graph); target++) {
            if (target == source || distances[source][target] == std::numeric_limits<double>::max())
                continue;
            Point_2 sourcePoint = graph[source];
            Point_2 targetPoint = graph[target];
            double stretch = distances[source][target] / euclidean_distance(sourcePoint, targetPoint);
            maxStretch = std::max(maxStretch, stretch);
        }
    }

    return maxStretch;
}

/* 
 * Function to compute the diameter of a t_spanner
 * (Algorithm by M.Farshi. PhD thesis: https://doi.org/10.6100/IR630219)
 */
int compute_t_spanner_diameter(const Graph &G, double t) {
    const size_t num_vertices = boost::num_vertices(G);

    std::vector<std::vector<double>> d_1(num_vertices, std::vector<double>(num_vertices, INF));
    std::vector<std::pair<Vertex, Vertex>> s_1;

    for (size_t i = 0; i < num_vertices; ++i) {
        d_1[i][i] = INF;
        for (size_t j = i + 1; j < num_vertices; ++j) {
            std::pair<Edge, bool> edge_pair = boost::edge(i, j, G);
            if (edge_pair.second) {
                d_1[i][j] = d_1[j][i] = boost::get(boost::edge_weight, G, edge_pair.first);
                s_1.emplace_back(i, j);
            } else {
                d_1[i][j] = d_1[j][i] = INF;
            }
        }
    }

    std::vector<std::vector<double>> d_k = d_1;
    std::vector<std::pair<Vertex, Vertex>> s_k = s_1;
    std::vector<std::pair<Vertex, Vertex>> s_next;
    std::vector<std::vector<double>> d_next(num_vertices, std::vector<double>(num_vertices, INF));

    int k = 1;
    bool flag = true;
    while (flag) {
        k++;
        d_next = d_k;
        flag = false;

        for (const auto &pair: s_k) {
            Vertex p = pair.first;
            Vertex q = pair.second;

            adjacency_iterator ai, ai_end;

            for (std::tie(ai, ai_end) = boost::adjacent_vertices(p, G); ai != ai_end; ++ai) {
                Vertex r = *ai;
                if (d_next[r][q] > d_k[p][q] + d_1[p][r]) {
                    d_next[r][q] = d_k[p][q] + d_1[p][r];
                    s_next.emplace_back(r, q);
                    if (d_next[r][q] > t * (euclidean_distance(G[r], G[q]))) {
                        flag = true;
                    }
                }
            }

            for (std::tie(ai, ai_end) = boost::adjacent_vertices(q, G); ai != ai_end; ++ai) {
                Vertex r = *ai;
                if (d_next[r][p] > d_k[p][q] + d_1[q][r]) {
                    d_next[r][p] = d_k[p][q] + d_1[q][r];
                    s_next.emplace_back(r, p);
                    if (d_next[r][p] > t * (euclidean_distance(G[r], G[p]))) {
                        flag = true;
                    }
                }
            }
        }

        d_k = d_next;
        s_k = s_next;
        s_next.clear();
    }
    return k;
}

/*
 * Constructing plane geometric spanner using the diametral method
 * https://doi.org/10.20382/jocg.v7i1a21
 */
std::vector<std::pair<Point_2, Point_2>> diametral_plane_spanner(const std::vector<Point_2> &points) {
    std::vector<Point_2> B = points;
    std::vector<std::pair<Point_2, Point_2>> E;

    for (std::size_t i = 0; i < B.size(); ++i) {
        std::size_t j = (i + 1) % B.size();
        E.emplace_back(B[i], B[j]);
    }

    while (B.size() >= 4) {
        Point_2 p;

        if (get_random_boolean_values())
            p = find_diametral_pair(B, true).first;
        else
            p = find_diametral_pair(B, true).second;

        Point_2 q, r;
        std::tie(q, r) = get_neighbors_on_convex_hull(B, p);

        E.emplace_back(q, r);

        auto it = std::find(B.begin(), B.end(), p);
        if (it != B.end())
            B.erase(it);
    }
    return E;
}

/*
 * Constructing plane geometric spanner using the approximate_diametral method
 * https://doi.org/10.20382/jocg.v7i1a21
 */
std::vector<std::pair<Point_2, Point_2>> approximate_diametral_plane_spanner(const std::vector<Point_2> &points) {
    std::vector<Point_2> B = points;
    std::vector<std::pair<Point_2, Point_2>> E;

    for (std::size_t i = 0; i < B.size(); ++i) {
        std::size_t j = (i + 1) % B.size();
        E.emplace_back(B[i], B[j]);
    }

    while (B.size() >= 4) {
        Point_2 p;

        if (get_random_boolean_values())
            p = find_approximate_diametral_pair(B, 112).first;
        else
            p = find_approximate_diametral_pair(B, 112).second;

        Point_2 q, r;
        std::tie(q, r) = get_neighbors_on_convex_hull(B, p);
        E.emplace_back(q, r);

        auto it = std::find(B.begin(), B.end(), p);
        if (it != B.end())
            B.erase(it);
    }
    return E;
}

// Constructing plane geometric spanner using the proposed method
std::pair<double, std::vector<std::pair<Point_2, Point_2>>> proposed_plane_spanner(const std::vector<Point_2> &points) {
    double myDilation = std::numeric_limits<double>::lowest();
    std::vector<Point_2> B = points;
    double perimeter = compute_perimeter(points);
    std::vector<std::pair<Point_2, Point_2>> E;

    for (std::size_t i = 0; i < B.size(); ++i) {
        std::size_t j = (i + 1) % B.size();
        E.emplace_back(B[i], B[j]);
    }

    while (B.size() >= 4) {
        Point_2 p;
        double t;
        std::tie(t, p) = find_proposed_t_good(B, perimeter);
        myDilation = std::max(myDilation, t);
        Point_2 q, r;
        std::tie(q, r) = get_neighbors_on_convex_hull(B, p);
        perimeter = perimeter - euclidean_distance(p, r) - euclidean_distance(p, q) + euclidean_distance(q, r);
        E.emplace_back(q, r);
        auto it = std::find(B.begin(), B.end(), p);
        if (it != B.end())
            B.erase(it);
    }
    return std::make_pair(myDilation, E);
}

/**
 * Function to store the graph in an adjacency list.
 *  A std::map is used to keep track of vertices to ensure that each point is added only once. This avoids adding duplicate vertices for the same point.
 *  Edges are added using the Euclidean distance as their weights.
 *  The distance matrix is initialized with std::numeric_limits<double>::infinity().
 *  Diagonal elements are set to 0.0 because the distance from a vertex to itself is zero.
 */
std::pair<Graph, double> build_graph(const std::vector<std::pair<Point_2, Point_2>> &spanner_edges) {
    Graph graph;

    // Map to keep track of unique vertices
    std::map<Point_2, Vertex> vertex_map;

    double weight = 0.0;

    // Add vertices and edges to the graph
    for (const auto &edge: spanner_edges) {
        Vertex v1, v2;

        // Check if the vertex already exists in the map
        if (vertex_map.find(edge.first) == vertex_map.end()) {
            v1 = add_vertex(edge.first, graph); // Add the vertex if it does not exist
            vertex_map[edge.first] = v1; // Map the point to the vertex descriptor
        } else {
            v1 = vertex_map[edge.first]; // Retrieve the existing vertex descriptor
        }

        if (vertex_map.find(edge.second) == vertex_map.end()) {
            v2 = add_vertex(edge.second, graph); // Add the vertex if it does not exist
            vertex_map[edge.second] = v2; // Map the point to the vertex descriptor
        } else {
            v2 = vertex_map[edge.second]; // Retrieve the existing vertex descriptor
        }

        double eu = euclidean_distance(edge.first, edge.second);
        weight += eu;
        // Add the edge with the weight being the Euclidean distance
        add_edge(v1, v2, eu, graph);
    }
    return std::make_pair(graph, weight);
}

// Function to find the maximum vertex degree and graph size.
std::pair<size_t, size_t> return_max_degree_and_size(const Graph &graph) {
    size_t numEdges = boost::num_edges(graph);
    size_t maxDegree = 0;

    for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
        Vertex v = *vp.first;
        size_t degree = out_degree(v, graph);
        if (degree > maxDegree) {
            maxDegree = degree; // Update the maximum degree found
        }
    }
    return std::make_pair(maxDegree, numEdges);
}

// Function to find the minimum spanning tree of input points
std::pair<double, double> compute_mst(const std::vector<Point_2> &points) {
    Graph g;
    std::vector<Vertex> vertices(points.size());

    for (size_t i = 0; i < points.size(); ++i) {
        vertices[i] = boost::add_vertex(points[i], g);
    }

    double t_distance = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            double distance = euclidean_distance(g[vertices[i]], g[vertices[j]]);
            boost::add_edge(vertices[i], vertices[j], distance, g);
            t_distance += distance;
        }
    }
    std::vector<Vertex> p(boost::num_vertices(g));
    prim_minimum_spanning_tree(g, &p[0]);

    double total_weight = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
        if (p[i] != i) {
            Edge e;
            bool exists;
            boost::tie(e, exists) = edge(p[i], i, g);
            if (exists) {
                total_weight += get(boost::edge_weight, g, e);
            }
        }
    }
    return std::make_pair(t_distance, total_weight);
}

// Function to read input points from the file.
std::vector<Point_2> read_points_from_file(const std::string &filename) {
    std::ifstream input(filename);
    std::vector<Point_2> points;
    double x, y;
    while (input >> x >> y) {
        points.emplace_back(x, y);
    }
    return points;
}

int main() {
    double diameter_entry_t = 1.88;
    //    *****  Read points from the file  *****
    std::string file_name = "1000";
    std::vector<Point_2> points = read_points_from_file("input_data/" + file_name + ".txt");

    //    *****  Compute MST of points  *****
    double c_g_weight, mst_weight;
    std::tie(c_g_weight, mst_weight) = compute_mst(points);
    std::cout << "File name: " << file_name << std::endl;
    std::cout << "Complete graph Weight: " << c_g_weight << "  MST Weight: " << mst_weight << std::endl;

    //    *****  Constructing a plane spanner with diametral approach  *****
    auto diametralStart = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<Point_2, Point_2>> diametral_spanner_edges = diametral_plane_spanner(points);
    auto diametralEnd = std::chrono::high_resolution_clock::now();

    // Calculate and print the duration
    std::chrono::duration<double> diametralDuration = diametralEnd - diametralStart;
    std::cout << "Diametral Plane Spanner Time : " << diametralDuration.count() << " seconds" << std::endl;

    std::string file_name1 = "output_spanner/ds_" + file_name + ".txt";
    write_spanner_to_file(diametral_spanner_edges, file_name1);

    //    *****  Constructing a plane spanner with approximate diametral approach  *****
    auto approximateDiametralStart = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<Point_2, Point_2>> approximate_diametral_spanner_edges = approximate_diametral_plane_spanner(points);
    auto approximateDiametralEnd = std::chrono::high_resolution_clock::now();

    // Calculate and print the duration
    std::chrono::duration<double> approximateDiametralDuration = approximateDiametralEnd - approximateDiametralStart;
    std::cout << "Approximate Diametral Plane Spanner Time : " << approximateDiametralDuration.count() << " seconds" << std::endl;

    std::string approximate_file_name1 = "output_spanner/ads_" + file_name + ".txt";
    write_spanner_to_file(approximate_diametral_spanner_edges, approximate_file_name1);

    //    *****  Constructing a plane spanner with proposed approach  *****
    double proposed_spanner_dilation;
    std::vector<std::pair<Point_2, Point_2>> proposed_spanner_edges;

    auto proposedStart = std::chrono::high_resolution_clock::now();
    std::tie(proposed_spanner_dilation, proposed_spanner_edges) = proposed_plane_spanner(points);
    auto proposedEnd = std::chrono::high_resolution_clock::now();

    // Calculate and print the duration
    std::chrono::duration<double> proposedDuration = proposedEnd - proposedStart;
    std::cout << "Proposed Plane Spanner Time : " << proposedDuration.count() << " seconds" << std::endl;
    std::cout << "Proposed dilation: " << proposed_spanner_dilation << std::endl;

    std::string file_name2 = "output_spanner/" + std::to_string(proposed_spanner_dilation) + "_ps_" + file_name + ".txt";
    write_spanner_to_file(proposed_spanner_edges, file_name2);

    //    *****  compute spanner weight, max degree and spanner size  *****
    double d_s_weight, p_s_weight, a_d_s_weight;
    Graph diametral_graph, approximate_graph, proposed_graph;

    std::tie(proposed_graph, p_s_weight) = build_graph(proposed_spanner_edges);
    std::tie(diametral_graph, d_s_weight) = build_graph(diametral_spanner_edges);
    std::tie(approximate_graph, a_d_s_weight) = build_graph(approximate_diametral_spanner_edges);

    std::cout << "Diametral spanner weight: " << d_s_weight << std::endl;
    std::cout << "Approximate Diametral spanner weight: " << a_d_s_weight << std::endl;
    std::cout << "Proposed spanner weight: " << p_s_weight << std::endl;

    int proposed_diameter = 0, diametral_diameter = 0, approximate_diameter = 0;
    size_t proposed_max_degree = 0, proposed_graph_size = 0, diametral_max_degree = 0, diametral_graph_size = 0, approximate_max_degree = 0, approximate_graph_size = 0;

    std::tie(proposed_max_degree, proposed_graph_size) = return_max_degree_and_size(proposed_graph);
    std::tie(diametral_max_degree, diametral_graph_size) = return_max_degree_and_size(diametral_graph);
    std::tie(approximate_max_degree, approximate_graph_size) = return_max_degree_and_size(approximate_graph);

    std::cout << "Proposed max degree: " << proposed_max_degree << std::endl;
    std::cout << "Proposed graph size: " << proposed_graph_size << std::endl;
    std::cout << "Diametral max degree: " << diametral_max_degree << std::endl;
    std::cout << "Diametral graph size: " << diametral_graph_size << std::endl;
    std::cout << "Approximate max degree: " << approximate_max_degree << std::endl;
    std::cout << "Approximate graph size: " << approximate_graph_size << std::endl;

    //    *****  compute stretch factor of Plane Spanners  *****
    auto diametral_stretch_factor_Start = std::chrono::high_resolution_clock::now();
    double diametral_stretch_factor = compute_stretch_factor(diametral_graph);
    auto diametral_stretch_factor_End = std::chrono::high_resolution_clock::now();

    auto approximate_stretch_factor_Start = std::chrono::high_resolution_clock::now();
    double approximate_stretch_factor = compute_stretch_factor(approximate_graph);
    auto approximate_stretch_factor_End = std::chrono::high_resolution_clock::now();

    auto proposed_stretch_factor_Start = std::chrono::high_resolution_clock::now();
    double proposed_stretch_factor = compute_stretch_factor(proposed_graph);
    auto proposed_stretch_factor_End = std::chrono::high_resolution_clock::now();

    // Calculate and print the duration
    std::chrono::duration<double> diametral_stretch_factor_Duration = diametral_stretch_factor_End - diametral_stretch_factor_Start;
    std::cout << "Stretch Factor Time of Diametral Plane Spanner: " << diametral_stretch_factor_Duration.count() << " seconds" << std::endl;

    std::chrono::duration<double> approximate_stretch_factor_Duration = approximate_stretch_factor_End - approximate_stretch_factor_Start;
    std::cout << "Stretch Factor Time of Approximate Diametral Plane Spanner: " << approximate_stretch_factor_Duration.count() << " seconds" << std::endl;

    std::chrono::duration<double> proposed_stretch_factor_Duration = proposed_stretch_factor_End - proposed_stretch_factor_Start;
    std::cout << "Stretch Factor Time of Proposed Plane Spanner: " << proposed_stretch_factor_Duration.count() << " seconds" << std::endl;

    std::cout << "Diametral Stretch Factor: " << diametral_stretch_factor << std::endl << "Approximate Diametral Stretch Factor: " << approximate_stretch_factor << std::endl << "Proposed Stretch Factor: " << proposed_stretch_factor << std::endl;

    //    *****  compute diameter of Plane Spanners  *****
    auto compute_t_spanner_diameter_Start = std::chrono::high_resolution_clock::now();
    proposed_diameter = compute_t_spanner_diameter(proposed_graph, diameter_entry_t);
    auto compute_t_spanner_diameter_End = std::chrono::high_resolution_clock::now();

    auto compute_t_spanner_diametral_diameter_Start = std::chrono::high_resolution_clock::now();
    diametral_diameter = compute_t_spanner_diameter(diametral_graph, diameter_entry_t);
    auto compute_t_spanner_diametral_diameter_End = std::chrono::high_resolution_clock::now();

    auto compute_approximate_diameter_Start = std::chrono::high_resolution_clock::now();
    approximate_diameter = compute_t_spanner_diameter(approximate_graph, diameter_entry_t);
    auto compute_approximate_diameter_End = std::chrono::high_resolution_clock::now();

    // Calculate and print the duration
    std::chrono::duration<double> compute_t_spanner_diameter_Duration = compute_t_spanner_diameter_End - compute_t_spanner_diameter_Start;
    std::cout << "compute proposed t spanner diameter time: " << compute_t_spanner_diameter_Duration.count() << " seconds" << std::endl;

    std::chrono::duration<double> compute_t_spanner_diametral_diameter_Duration = compute_t_spanner_diametral_diameter_End - compute_t_spanner_diametral_diameter_Start;
    std::cout << "compute diametral t spanner diameter time: " << compute_t_spanner_diametral_diameter_Duration.count()<< " seconds" << std::endl;

    std::chrono::duration<double> compute_approximate_diameter_Duration = compute_approximate_diameter_End - compute_approximate_diameter_Start;
    std::cout << "compute approximate diametral t spanner diameter time: " << compute_approximate_diameter_Duration.count() << " seconds" << std::endl;

    std::cout << "Proposed diameter: " << proposed_diameter << std::endl;
    std::cout << "Diametral diameter: " << diametral_diameter << std::endl;
    std::cout << "Approximate diameter: " << approximate_diameter << std::endl;
    std::cout << "The END" << std::endl;

    return 0;
}

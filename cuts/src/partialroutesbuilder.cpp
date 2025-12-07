#include "partialroutesbuilder.h"
#include "graphviewer.h"

PartialRoutesBuilder::PartialRoutesBuilder(const SVRPInstance &instance,
                                           const Params &params,
                                           const EdgeValueMap &xValue,
                                           const NodeValueMap &recourseValue)
    : instance(instance), params(params), xValue(xValue),
      recourseValue(recourseValue), supportXValue(supportGraph),
      edgeDepotValue(supportGraph), treeVerticesMap(tree),
      treeCutVertexMap(tree), treeRecourseCostMap(tree), treeNodeXValue(tree),
      treeDemand(tree) {
    buildPartialRoutes();
}

void PartialRoutesBuilder::buildPartialRoutes() {
    // Construct support graph.
    buildSupportGraph();

    // Check if support graph is empty or has only one node.
    int numberOfNodes = 0;
    for (NodeIt v(supportGraph); v != INVALID; ++v) {
        numberOfNodes++;

        if (numberOfNodes >= 2) {
            break;
        }
    }

    if (numberOfNodes <= 1) {
        return;
    }

    // Build block-cut trees from the support graph.
    buildBlockCutTree();

    // Build the partial routes.
    // First collect all leaves and singletons.
    std::vector<Node> leaves;
    std::vector<Node> singletons;

    for (NodeIt v(tree); v != INVALID; ++v) {
        int numberOfIncEdges = 0;
        for (IncEdgeIt e(tree, v); e != INVALID; ++e) {
            numberOfIncEdges++;
        }

        if (numberOfIncEdges == 1) {
            leaves.push_back(v);
        }

        if (numberOfIncEdges == 0) {
            singletons.push_back(v);
        }
    }

    // Add partial route for singletons.
    for (Node v : singletons) {
        PartialRoute partialRoute;
        pushNodeIntoPartialRoute(partialRoute, v);
        partialRoute.isRoute = isRoute(partialRoute);
        checkAndSortPartialRoute(partialRoute);

        // check if 1 + (x(H) - |H| + 1) >= eps
        assert(static_cast<int>(partialRoute.entries.size()) >= 1);
        if (partialRoute.entries[0].xInside -
                static_cast<int>(partialRoute.entries[0].vertices.size()) >=
            -2.0 + EpsForIntegrality) {
            partialRoutes.push_back(partialRoute);
        }
    }

    // For each leaf, call the dfs.
    // (The cast is due to the case leaves.size() = 0.)
    for (int i = 0; i < static_cast<int>(leaves.size()) - 1; i++) {
        DFS dfs(tree);
        dfs.run(leaves[i]);

        for (size_t j = i + 1; j < leaves.size(); j++) {
            if (dfs.reached(leaves[j])) {
                // Contruct path.
                std::vector<Node> path = {leaves[j]};
                Node prev = dfs.predNode(leaves[j]);
                while (prev != INVALID) {
                    path.push_back(prev);
                    prev = dfs.predNode(prev);
                }

                // Add partial route corresponding to the path from leaf i
                // to leaf j.
                PartialRoute partialRoute;
                double incrementPartialRoute = 0.0;
                double totalDemand = 0.0;
                std::vector<int> vertices;
                for (Node &treeNode : path) {
                    totalDemand += treeDemand[treeNode];
                    pushNodeIntoPartialRoute(partialRoute, treeNode);
                    if (!treeVerticesMap[treeNode].empty()) {
                        incrementPartialRoute +=
                            partialRoute.entries.back().xInside +
                            partialRoute.entries.back().xLeft;
                    }
                    for (const int id : treeVerticesMap[treeNode]) {
                        vertices.push_back(id);
                    }
                }

                partialRoute.isRoute = isRoute(partialRoute);
                checkAndSortPartialRoute(partialRoute);

                int kH = static_cast<int>(
                    std::ceil(totalDemand / instance.capacity));

                // By now, `incrementPartialRoute` has x(H), so here we
                // check if 1 + (x(H) - |H| + \lceil d(S) / C \rceil) >= eps.
                if (incrementPartialRoute + kH -
                            static_cast<int>(vertices.size()) >=
                        -1.0 + EpsForIntegrality &&
                    customerSets.find(vertices) != customerSets.end()) {
                    customerSets.insert(vertices);
                }

                // Similarly, here we check if 1 + (x(H) - |H| + 1) >= eps.
                if (incrementPartialRoute - static_cast<int>(vertices.size()) >=
                    -2.0 + EpsForIntegrality) {
                    partialRoutes.push_back(partialRoute);
                }
            }
        }
    }
}

void PartialRoutesBuilder::buildSupportGraph() {
    // Add nodes.
    for (int i = 0; i < instance.n; i++) {
        supportGraph.addNode();
    }

    // Add edges in the support.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        Node u2 = supportGraph.nodeFromId(instance.g.id(u));
        Node v2 = supportGraph.nodeFromId(instance.g.id(v));
        assert(u2 != INVALID && v2 != INVALID);

        if (xValue[e] >= EpsForIntegrality) {
            Edge e2 = supportGraph.addEdge(u2, v2);
            supportXValue[e2] += xValue[e];
        }
    }

    // Get value of depot edges.
    Node depot = instance.g.nodeFromId(0);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (v == depot) {
            continue;
        }

        Node v2 = supportGraph.nodeFromId(instance.g.id(v));
        Edge depotEdge = findEdge(instance.g, v, depot);
        assert(v2 != INVALID && depotEdge != INVALID);
        edgeDepotValue[v2] = xValue[depotEdge];
    }

    // Find connected components of the support graph.
    NodeIntMap componentsMap(supportGraph);
    int numberOfComponents = connectedComponents(supportGraph, componentsMap);

    // Count the sum of the "edge depot values" for each component.
    // Also get the nodes and total demand in each component.
    std::vector<float> componentEdgeDepotValues(numberOfComponents, 0.0);
    std::vector<std::vector<int>> componentVertices(numberOfComponents,
                                                    std::vector<int>());
    std::vector<double> componentDemand(numberOfComponents, 0.0);
    for (NodeIt v(supportGraph); v != INVALID; ++v) {
        int componentId = componentsMap[v];
        assert(componentId >= 0 && componentId < numberOfComponents);
        componentVertices[componentId].push_back(supportGraph.id(v));
        componentDemand[componentId] +=
            instance.demand[instance.g.nodeFromId(supportGraph.id(v))];
        componentEdgeDepotValues[componentsMap[v]] += edgeDepotValue[v];
    }

    // Add customer sets for each component.
    for (int i = 0; i < numberOfComponents; i++) {
        if (static_cast<int>(componentVertices[i].size()) == 0 ||
            componentDemand[i] <= EpsForIntegrality) {
            continue;
        }

        customerSets.insert(componentVertices[i]);
    }

    // Delete components whose total edge depot value is different than two.
    Node newDepot = supportGraph.nodeFromId(instance.depot);
    std::vector<Node> toErase = {newDepot};
    if (!params.scenarioOptimalPRCuts) {
        for (NodeIt v(supportGraph); v != INVALID; ++v) {
            if (v == newDepot) {
                continue;
            }

            assert(componentsMap[v] >= 0 &&
                   componentsMap[v] < numberOfComponents);
            if (std::abs(componentEdgeDepotValues[componentsMap[v]] - 2.0) >=
                EpsForIntegrality) {
                toErase.push_back(v);
            }
        }
    }

    for (Node v : toErase) {
        supportGraph.erase(v);
    }
}

void PartialRoutesBuilder::buildBlockCutTree() {
    // Find the biconnected components.
    EdgeIntMap biconnectedEdgeMap(supportGraph);
    int numberOfBiComponents =
        biNodeConnectedComponents(supportGraph, biconnectedEdgeMap);

    // Find the cut nodes.
    NodeBoolMap boolCutNodesMap(supportGraph);
    biNodeConnectedCutNodes(supportGraph, boolCutNodesMap);

    // Attribute ids to nodes in biconnected components and to cut nodes.
    // Ids of bicomponents start at zero and id of cut nodes start at
    // `numberOfBiComponents`.
    int numberOfCutNodes = 0;
    NodeIntMap nodeIdMap(supportGraph);
    for (NodeIt v(supportGraph); v != INVALID; ++v) {
        if (boolCutNodesMap[v] ||
            edgeDepotValue[v] >= 1.0 - EpsForIntegrality) {
            nodeIdMap[v] = numberOfBiComponents + numberOfCutNodes;
            numberOfCutNodes++;
        } else {
            IncEdgeIt e(supportGraph, v);
            if (e == INVALID) {
                // Single isolated node.
                assert(edgeDepotValue[v] >= 2.0 - EpsForIntegrality);
                nodeIdMap[v] = numberOfBiComponents + numberOfCutNodes;
                numberOfCutNodes++;
            } else {
                nodeIdMap[v] = biconnectedEdgeMap[e];
            }
        }
    }

    // Add nodes to the tree.
    for (int i = 0; i < numberOfBiComponents + numberOfCutNodes; i++) {
        Node v = tree.addNode();
        treeRecourseCostMap[v] = INFINITY;
        treeNodeXValue[v] = 0.0;
        treeDemand[v] = 0.0;

        if (tree.id(v) < numberOfBiComponents) {
            treeCutVertexMap[v] = false;
        } else {
            treeCutVertexMap[v] = true;
        }
    }

    // Add edges to the tree.
    for (EdgeIt e(supportGraph); e != INVALID; ++e) {
        Node component = tree.nodeFromId(biconnectedEdgeMap[e]);
        int nodeId1 = nodeIdMap[supportGraph.u(e)];
        int nodeId2 = nodeIdMap[supportGraph.v(e)];
        Node node1 = tree.nodeFromId(nodeId1);
        Node node2 = tree.nodeFromId(nodeId2);

        if (nodeId1 >= numberOfBiComponents &&
            findEdge(tree, component, node1) == INVALID) {
            tree.addEdge(component, node1);
        }
        if (nodeId2 >= numberOfBiComponents &&
            findEdge(tree, component, node2) == INVALID) {
            tree.addEdge(component, node2);
        }
    }

    // Update tree information.
    for (NodeIt node(supportGraph); node != INVALID; ++node) {
        int id = supportGraph.id(node);
        Node originalNode = instance.g.nodeFromId(id);
        int componentId =
            nodeIdMap[node]; // can be either cut node id or component id.

        Node treeNode = tree.nodeFromId(componentId);
        treeVerticesMap[treeNode].push_back(id);
        treeDemand[treeNode] += instance.demand[originalNode];
        treeRecourseCostMap[treeNode] =
            std::min(treeRecourseCostMap[treeNode],
                     instance.getEdgeRecourseCost(originalNode));
    }

    for (NodeIt treeNode(tree); treeNode != INVALID; ++treeNode) {
        double nodeXValue = 0.0;
        for (int i = 0; i < static_cast<int>(treeVerticesMap[treeNode].size());
             i++) {
            for (int j = i + 1;
                 j < static_cast<int>(treeVerticesMap[treeNode].size()); j++) {
                Node u = instance.g.nodeFromId(treeVerticesMap[treeNode][i]);
                Node v = instance.g.nodeFromId(treeVerticesMap[treeNode][j]);
                Edge e = findEdge(instance.g, u, v);
                assert(u != INVALID && v != INVALID && e != INVALID);
                nodeXValue += xValue[e];
            }
        }
        treeNodeXValue[treeNode] = nodeXValue;
    }
}

// This functions push a set S into a partial route H = (S_1, \ldots,
// S_\ell) and returns x(S_\ell, S) + x(S) - |S|.
void PartialRoutesBuilder::pushNodeIntoPartialRoute(PartialRoute &partialRoute,
                                                    Node treeNode) {
    if (treeVerticesMap[treeNode].empty()) {
        return;
    }

    double xInside = treeNodeXValue[treeNode];
    double xRight = 0.0;
    if (!partialRoute.entries.empty()) {
        for (const int id1 : partialRoute.entries.back().vertices) {
            for (const int id2 : treeVerticesMap[treeNode]) {
                Node u = instance.g.nodeFromId(id1);
                Node v = instance.g.nodeFromId(id2);
                Edge e = findEdge(instance.g, u, v);
                xRight += xValue[e];
            }
        }

        assert(partialRoute.entries.back().xRight == 0.0);
        partialRoute.entries.back().xRight = xRight;
    }

    partialRoute.entries.push_back({
        treeVerticesMap[treeNode],     // vertices
        treeRecourseCostMap[treeNode], // recourseCost
        xInside,                       // xInside
        xRight,                        // xLeft
        0.0                            // xRight
    });
}

bool PartialRoutesBuilder::isRoute(const PartialRoute &partialRoute) {
    for (const PartialRouteEntry &entry : partialRoute.entries) {
        if (static_cast<int>(entry.vertices.size()) > 1) {
            return false;
        }
    }

    return true;
}

// This checks if the partial route is ok and sort the nodes in an unstructured
// set by their recourse costs.
void PartialRoutesBuilder::checkAndSortPartialRoute(
    PartialRoute &partialRoute) {
    std::unordered_set<int> visited;

    for (size_t i = 0; i < partialRoute.entries.size(); i++) {
        if (partialRoute.entries[i].vertices.size() > 1) {
            // Check partial route condition.
            assert(i == static_cast<int>(partialRoute.entries.size()) - 1 ||
                   static_cast<int>(
                       partialRoute.entries[i + 1].vertices.size()) == 1);

            // Sort vertices by their recourse cost.
            std::sort(partialRoute.entries[i].vertices.begin(),
                      partialRoute.entries[i].vertices.end(),
                      [&](int id1, const int id2) {
                          Node u = instance.g.nodeFromId(id1);
                          Node v = instance.g.nodeFromId(id2);
                          assert(u != INVALID && v != INVALID);
                          return instance.getEdgeRecourseCost(u) <
                                 instance.getEdgeRecourseCost(v);
                      });
        }

        // Recourse should be minimum and no vertex should be repeated.
        double minRecourse = DBL_MAX;
        for (int id : partialRoute.entries[i].vertices) {
            assert(visited.count(id) == 0);
            Node v = instance.g.nodeFromId(id);
            minRecourse =
                std::min(minRecourse, instance.getEdgeRecourseCost(v));
            visited.insert(id);
        }
        assert(std::abs(minRecourse - partialRoute.entries[i].recourseCost) <=
               1e-6);
    }
}
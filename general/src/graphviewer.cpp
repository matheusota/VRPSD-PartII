#include "graphviewer.h"

// just some colors that we can use
namespace graphviewer {
std::vector<std::string> colors = {
    "Gray",   "Red",  "Blue",  "Orange", "Green", "Pink", "Brown", "Purple",
    "Yellow", "Cyan", "Olive", "Beige",  "Lime",  "Teal", "Navy",  "Lavender"};

void viewVRPSolution(const SVRPInstance &instance, const EdgeValueMap &x,
                     const std::string title) {
    // assign names to the nodes
    NodeStringMap nodeNames(instance.g);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        nodeNames[v] = "\"" + to_string(instance.g.id(v)) + "(" +
                       to_string(static_cast<int>(instance.demand[v])) + ")" +
                       "\"";
    }

    // set graph attributes for the visualizer
    GraphAttributes GA(instance.g, nodeNames, instance.posx, instance.posy);
    GA.SetLabel(title);

    GA.SetDefaultNodeAttrib("fixedsize = true width = 1 height = 1 style = "
                            "filled fixedsize = true fontsize = 20");
    GA.SetDefaultEdgeAttrib("fontsize = 20");

    //    GA.SetDefaultNodeAttrib("color = Gray shape = ellipse fontsize = 35");
    //    GA.SetDefaultNodeAttrib("fixedsize = true color = Gray shape = circle
    //    style = filled fontsize = 0"); GA.SetDefaultEdgeAttrib("color =
    //    Gray");

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        GA.SetColor(v, colors[0]);
    }

    // color arcs according to the solution
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (x[e] > 0.0001) {
            GA.SetColor(e, colors[1]);

            stringstream stream;
            stream << fixed << setprecision(2) << x[e];
            string s = stream.str();

            GA.SetAttrib(e, "label = " + s);
        } else {
            GA.SetColor(e, "Invis");
        }
    }

    // view graph
    GA.View();
}

void viewVRPSolution(const SVRPInstance &instance, const SVRPSolution &solution,
                     const std::string title) {

    // create digraph for visualization
    Digraph g;
    DNodeValueMap posx(g);
    DNodeValueMap posy(g);

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        g.addNode();
    }

    // assign names and positions to the nodes
    DNodeStringMap nodeNames(g);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        DNode v2 = g.nodeFromId(instance.g.id(v));
        nodeNames[v2] = "\"" + to_string(instance.g.id(v)) + "(" +
                        to_string(static_cast<int>(instance.demand[v])) + ")" +
                        "\"";
        posx[v2] = instance.posx[v];
        posy[v2] = instance.posy[v];
    }

    // set graph attributes for the visualizer
    DigraphAttributes GA(g, nodeNames, posx, posy);
    GA.SetLabel(title);

    GA.SetDefaultDNodeAttrib("fixedsize = true width = 2 height = 2 style = "
                             "filled fixedsize = true fontsize = 30");
    GA.SetDefaultArcAttrib("penwidth = 4");

    //    GA.SetDefaultNodeAttrib("color = Gray shape = ellipse fontsize = 35");
    //    GA.SetDefaultNodeAttrib("fixedsize = true color = Gray shape = circle
    //    style = filled fontsize = 0"); GA.SetDefaultEdgeAttrib("color =
    //    Gray");

    for (DNodeIt v(g); v != INVALID; ++v) {
        GA.SetColor(v, colors[0]);
    }

    DNode depot = g.nodeFromId(instance.depot);
    Arc e;
    for (size_t i = 0; i < solution.routes.size(); i++) {
        e = g.addArc(depot, g.nodeFromId(solution.routes[i][0]));
        GA.SetColor(e, colors[i + 1]);
        for (size_t j = 1; j < solution.routes[i].size(); j++) {
            e = g.addArc(g.nodeFromId(solution.routes[i][j - 1]),
                         g.nodeFromId(solution.routes[i][j]));
            GA.SetColor(e, colors[i + 1]);
        }
        e = g.addArc(
            g.nodeFromId(solution.routes[i][solution.routes[i].size() - 1]),
            depot);
        GA.SetColor(e, colors[i + 1]);
    }

    // view graph
    GA.View();
}

void viewSupportWithNodeAndEdgeMap(const Graph &g, const EdgeValueMap &x,
                                   const NodePosMap &posx,
                                   const NodePosMap &posy,
                                   const NodeBoolMap &nodeMap,
                                   const EdgeIntMap &edgeMap) {

    // Assign names to the nodes.
    NodeStringMap nodeNames(g);
    for (NodeIt v(g); v != INVALID; ++v) {
        nodeNames[v] = "\"" + to_string(g.id(v)) + "\"";
    }

    // set graph attributes for the visualizer
    GraphAttributes GA(g, nodeNames, posx, posy);

    GA.SetDefaultNodeAttrib("fixedsize = true width = 1 height = 1 style = "
                            "filled fixedsize = true fontsize = 20");
    GA.SetDefaultEdgeAttrib("fontsize = 20");

    //    GA.SetDefaultNodeAttrib("color = Gray shape = ellipse fontsize = 35");
    //    GA.SetDefaultNodeAttrib("fixedsize = true color = Gray shape = circle
    //    style = filled fontsize = 0"); GA.SetDefaultEdgeAttrib("color =
    //    Gray");

    for (NodeIt v(g); v != INVALID; ++v) {
        if (nodeMap[v]) {
            GA.SetColor(v, colors[1]);
        } else {
            GA.SetColor(v, colors[0]);
        }
    }

    // color arcs according to the solution
    for (EdgeIt e(g); e != INVALID; ++e) {
        if (x[e] > 0.0001) {
            GA.SetColor(e, colors[edgeMap[e] % 16]);

            stringstream stream;
            stream << fixed << setprecision(2) << x[e];
            string s = stream.str();

            GA.SetAttrib(e, "label = " + s);
        } else {
            GA.SetColor(e, "Invis");
        }
    }

    // view graph
    GA.View();
}
} // namespace graphviewer

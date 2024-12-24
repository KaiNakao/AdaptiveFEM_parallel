#pragma once
#include <vector>
#include <iostream>
#include <set>

struct Edge {
    std::set<int> nodes;

    Edge() {
        nodes = {};
    }
    
    bool operator<(const Edge &rhs) const {
        return nodes < rhs.nodes;
    }
};


struct Face {
    std::set<int> nodes;
    Edge marked_edge;

    Face() {
        nodes = {};
        marked_edge = {};
    }
};


struct Tetrahedron {
    std::set<int> nodes;
    std::vector<Face> faces;
    Edge refinement_edge;
    bool flag;
    int material;
    bool original;

    Tetrahedron() {
        nodes = {};
        faces = {};
        refinement_edge = {};
        flag = false;
        material = -1;
        original = false;
    }

    bool operator<(const Tetrahedron &rhs) const {
        return nodes < rhs.nodes;
    }

    inline void print() const {
        std::cout << "tetrahedron nodes: ";
        for (int node : nodes) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
        std::cout << "faces: " << std::endl;
        for (const auto &face : faces) {
            std::cout << " face nodes: ";
            for (int node : face.nodes) {
                std::cout << node << " ";
            }
            std::cout << std::endl;
            std::cout << " marked edge nodes: ";
            for (int node : face.marked_edge.nodes) {
                std::cout << node << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "refinement edge nodes: ";
        for (int node : refinement_edge.nodes) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    inline void check() const {
        for (const auto &face : faces) {
            if (face.nodes.size() != 3) {
                std::cout << "face nodes size is not 3" << std::endl;
                std::exit(1);
            }
            if (face.marked_edge.nodes.size() != 2) {
                std::cout << "marked edge nodes size is not 2" << std::endl;
                std::exit(1);
            }
            for (int node_id : face.nodes) {
                if (nodes.count(node_id) == 0) {
                    std::cout << "face node not in tetrahedron nodes" << std::endl;
                    std::exit(1);
                }
            }
            for (int node_id : face.marked_edge.nodes) {
                if (face.nodes.count(node_id) == 0) {
                    std::cout << "marked edge node not in face nodes" << std::endl;
                    std::exit(1);
                }
            }
        }
        for (int node_id : refinement_edge.nodes) {
            if (nodes.count(node_id) == 0) {
                std::cout << "refinement edge node not in tetrahedron nodes" << std::endl;
                std::exit(1);
            }
        }
    }
};
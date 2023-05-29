use itertools::izip;
use nauty_pet::graph::CanonGraph;
use petgraph::Undirected;
use std::collections::hash_map;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::fs::File;
use std::io::{Write};


fn import_lattice() -> (neighbors: HashMap<usize, Vec<usize>>, )

fn vsimple(
    vertices: &HashMap<usize, Vec<usize>>,
    cluster: &mut Vec<usize>,
    neighbors: &mut Vec<usize>,
    guarding_set: &HashSet<usize>,
    size: usize,
    all_clusters: &mut Vec<Vec<usize>>,
) -> bool {
    if cluster.len() == size {
        all_clusters.push(cluster.clone());
        return true;
    };

    let mut has_int_leaf = false;

    let mut new_guarding_set = guarding_set.clone();

    while !neighbors.is_empty() {
        let neighbor = neighbors.pop().unwrap();

        if vertices.contains_key(&neighbor) {
            cluster.push(neighbor);

            let mut new_neighbors = neighbors.clone();

            for vertex in vertices[&neighbor].iter() {
                if !cluster.contains(vertex)
                    & !new_guarding_set.contains(vertex)
                    & !new_neighbors.contains(vertex)
                {
                    new_neighbors.push(*vertex);
                };
            }

            if vsimple(
                vertices,
                cluster,
                &mut new_neighbors,
                &new_guarding_set,
                size,
                all_clusters,
            ) {
                cluster.pop();
                has_int_leaf = true;
            } else {
                cluster.pop();
                return has_int_leaf;
            };

            new_guarding_set.insert(neighbor);

            if (vertices.len() - new_guarding_set.len()) < size {
                return has_int_leaf;
            };
        };
    }

    has_int_leaf
}

///
/// Wrapper for the vsimple algorithm, takes in the hashmap that
/// links each vertex to its edges, the final size of the clusters
/// that are under consideration, the starting vertex (usually the
/// center vertex, but in cases of breaking apart the entire cluster
/// it will be all the vertices) and a vector of all the clusters
///
fn enumerate(
    vertices: &HashMap<usize, Vec<usize>>, // Hash map of the form {Vertex: [Neighbors], ...}
    size: usize, // Maximum size of the cluster
    starting_vertices: &Vec<usize>, // Vector of the form [Vertex0, Vertex1, ...]
    mut all_clusters: &mut Vec<Vec<usize>>) // Vector of the form [[Vertex0, Vertex1, ...], [Vertex0, Vertex2, ...], ...]
    {
    // Declare an empty guarding set
    let mut guarding_set = HashSet::<usize>::new();

    // Start looping over all the possible starting vertices
    for vertex in starting_vertices {
        // Find the neighbors of the vertex from the input hashmap
        // then filter out all the neighbors that are in the guarding set
        let mut neighbors = vertices[vertex]
            .iter()
            .cloned()
            .filter(|neighbor| !guarding_set.contains(neighbor))
            .collect();

        // Seed the initial cluster with just the starting vertex
        let mut starting_cluster = vec![*vertex];
        // Call vsimple, it will edit all_clusters and once it is done, we will be able to
        // look at the clusters that are generated simulataneously
        vsimple(
            vertices,
            &mut starting_cluster,
            &mut neighbors,
            &guarding_set,
            size,
            &mut all_clusters,
        );
        // Add the starting vertex into the guarding set
        guarding_set.insert(*vertex);
    }
}

fn main(){
    use std::env;
    let args: Vec<_> = env::args().collect();

    // Import Data as a JSON file here instead of generating it here
    let nlce_type: String = args[1].parse().unwrap();
    let nlce_directory = format!("./NLCE_Data/{}", nlce_type);

    let mut all_clusters = Vec::new();

    enumerate(&cluster_set, cluster_size, &vec![start], &mut all_clusters);

}

#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]

use nauty_pet::graph::CanonGraph;
use petgraph::Undirected;
use std::hash::{Hash, Hasher};
use std::collections::{HashSet, HashMap};
use std::collections::hash_map;
use itertools::izip;

fn _print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn cluster_to_iso_sym(cluster: &[usize], edge_list: &[Vec<usize>], iso_types: &[Vec<u8>], sym_types: &[Vec<u32>]) -> (Vec<(usize, usize, u8)>, Vec<usize>){
    // Takes cluster (sorted) and returns the corresponding ismorphic edge list and vertex type list
    let mut uncanon_list: Vec<(usize, usize, u8)> = vec![];
    let mut vertex_type: Vec<usize> = vec![];
    for (index, &vertex) in cluster.iter().enumerate(){
        vertex_type.push(0);
        for (&edge, &weight, &direction) in izip!(&edge_list[vertex], &iso_types[vertex], &sym_types[vertex]){
            let edge_position = cluster.iter().position(|&x| x == edge);
            match edge_position {
                Some(x) => {
                    vertex_type[index] += 2usize.pow(direction);
                    if x > index {
                        uncanon_list.push((index, x, weight));
                    }
                },
                _ => continue,
            }
        }
    }
    (uncanon_list, vertex_type)
}

fn cluster_to_iso(cluster: &[usize], edge_list: &[Vec<usize>], iso_types: &[Vec<u8>]) -> Vec<(usize, usize, u8)>{
    // Takes cluster (unsorted) and returns the corresponding isomorphic edge list
    let mut uncanon_list: Vec<(usize, usize, u8)> = vec![];
    for (index, &vertex) in cluster.iter().enumerate(){
        for (&edge, &weight) in edge_list[vertex].iter().zip(iso_types[vertex].iter()){
            let edge_position = cluster.iter().position(|&x| x == edge);
            match edge_position {
                Some(x) =>
                    if x > index {
                        uncanon_list.push((index, x, weight));
                    },
                _ => continue,
            }
        }
    }
    uncanon_list
}

fn iso_to_hash(uncanon_list: &[(usize, usize, u8)]) -> usize {
    let canon_graph = CanonGraph::<(), u8, Undirected, usize>::from_edges(uncanon_list);
    let mut graph_hasher = hash_map::DefaultHasher::new();
    canon_graph.hash(&mut graph_hasher);
    graph_hasher.finish() as usize
}

fn sym_to_hash(vertex_type_cluster: &[usize]) -> usize {
    let mut sym_hasher = hash_map::DefaultHasher::new();
    vertex_type_cluster.hash(&mut sym_hasher);
    sym_hasher.finish() as usize
}

fn add_cluster(cluster: HashSet<usize>,
               edge_list: &[Vec<usize>],
               iso_types: &[Vec<u8>],
               sym_types: &[Vec<u32>],
               graph_multiplicity: &mut HashMap<usize, usize>,
               subgraph_multiplicity: &mut HashMap<usize, HashMap<usize, usize>>,
               graph_bond_info: &mut HashMap<usize, Vec<(usize, usize, u8)>>,
               sym_hash_set: &mut HashSet<usize>) {

    let mut cluster_vec: Vec<usize> = cluster.into_iter().collect();
    cluster_vec.sort();

    let (cluster_iso_list, vertex_type_cluster) = cluster_to_iso_sym(&cluster_vec, edge_list, iso_types, sym_types);

    let iso_hash = iso_to_hash(&cluster_iso_list);
    let sym_hash = sym_to_hash(&vertex_type_cluster);

    if !graph_multiplicity.contains_key(&iso_hash) {

        sym_hash_set.insert(sym_hash);
        graph_multiplicity.insert(iso_hash, 1);
        graph_bond_info.insert(iso_hash, cluster_iso_list);

        //let mut subgraph_func = |subcluster: HashSet<usize>| add_subcluster(subcluster, edge_list, iso_types, iso_hash, subgraph_multiplicity);

        //for size in 2..cluster_vec.len() {
        //    enumerate(&restricted_set, size, &cluster_vec, &mut subgraph_func);
        //}

    } else if !sym_hash_set.contains(&sym_hash) {
        sym_hash_set.insert(sym_hash);
        graph_multiplicity.insert(iso_hash, graph_multiplicity[&iso_hash] + 1);
    };

}

fn add_subcluster(cluster: HashSet<usize>, edge_list: &[Vec<usize>], iso_types: &[Vec<u8>], iso_hash: usize, subgraph_multiplicity: &mut HashMap<usize, HashMap<usize, usize>>) {
    let cluster_vec: Vec<usize> = cluster.into_iter().collect();

    let subcluster_iso_list = cluster_to_iso(&cluster_vec, edge_list, iso_types);
    let sub_iso_hash = iso_to_hash(&subcluster_iso_list);
    subgraph_multiplicity.entry(iso_hash).and_modify(|subgraph_info| {subgraph_info.entry(sub_iso_hash).and_modify(|counter| *counter += 1).or_insert(1);})
                                         .or_insert(HashMap::from([(sub_iso_hash, 1)]));

}

fn vsimple_vec(edge_list: &[Vec<usize>],
           subgraph: &mut HashSet<usize>,
           neighbors: &mut Vec<usize>,
           guarding_set: &HashSet<usize>,
           size: usize,
           graph_func: &mut dyn FnMut(HashSet<usize>)) ->  bool {

    if subgraph.len() == size {
        println!("{:?}",subgraph);
        graph_func(subgraph.clone());
        return true;
    };

    let mut has_int_leaf = false;

    let mut guarding_set_clone = guarding_set.clone();

        println!("{:?}", neighbors);
    while !neighbors.is_empty() {
        let neighbor = neighbors.pop().unwrap();
        subgraph.insert(neighbor);

        let mut new_neighbors = neighbors.clone();
        let full_guard: HashSet<usize> = subgraph.union(&guarding_set_clone).cloned().collect();
        let add_neighbors = edge_list[neighbor].iter().cloned().filter(|&vertex| !full_guard.contains(&vertex) & !neighbors.contains(&vertex));

            new_neighbors.extend(add_neighbors);

            if vsimple_vec(edge_list, subgraph, &mut new_neighbors, &guarding_set_clone, size, graph_func) {
                subgraph.remove(&neighbor);
                has_int_leaf = true;
            } else{
                subgraph.remove(&neighbor);
                return has_int_leaf;
            };

            guarding_set_clone.insert(neighbor);

            if (edge_list.len() - guarding_set_clone.len()) < size {
                return has_int_leaf;
            };
    };

    has_int_leaf
}

fn enumerate_vec(edge_list: &[Vec<usize>], size: usize, starting_vertices: &Vec<usize>, mut graph_func: &mut dyn FnMut(HashSet<usize>)) {

    let mut guarding_set = HashSet::<usize>::new();

    for vertex in starting_vertices {
        println!("hjaepd");
        let mut neighbors = edge_list[*vertex].iter().cloned().filter(|vertex| !guarding_set.contains(vertex)).collect();
        let mut starting_subgraph = HashSet::from([*vertex]);
        vsimple_vec(edge_list, &mut starting_subgraph, &mut neighbors, &guarding_set, size, &mut graph_func);
        guarding_set.insert(*vertex);
    };
}

fn main() {
    let cluster = vec![0, 1, 2, 3];
    let cluster_set: HashMap<usize, HashSet<usize>> = HashMap::from([(0, HashSet::from([1, 2, 3])),
                                                                     (1, HashSet::from([0, 3, 2])),
                                                                     (2, HashSet::from([3, 1, 0])),
                                                                     (3, HashSet::from([0, 1, 2]))]);
    let edge_list: Vec<Vec<usize>> = vec![vec![1, 2, 3], vec![0, 3, 2], vec![3, 1, 0], vec![0, 1, 2]];
    let iso_types: Vec<Vec<u8>> = vec![vec![1, 1, 2], vec![1, 1, 2], vec![1, 2, 1], vec![2, 1, 1]];
    let sym_types: Vec<Vec<u32>> = vec![vec![0, 2, 1], vec![4, 2, 3], vec![0, 7, 6], vec![5, 6, 4]];

    let mut graph_mult = HashMap::<usize, usize>::new();
    let mut subgraph_mult = HashMap::<usize, HashMap<usize, usize>>::new();
    let mut graph_bond = HashMap::<usize, Vec<(usize, usize, u8)>>::new();
    let mut sym_hash = HashSet::<usize>::new();

    let mut graph_func = |cluster: HashSet<usize>| add_cluster(cluster, &edge_list, &iso_types, &sym_types, &mut graph_mult, &mut subgraph_mult, &mut graph_bond, &mut sym_hash);

    enumerate_vec(&edge_list, 3, &cluster, &mut graph_func);
    println!("{:?}", graph_mult);
    println!("{:?}", graph_bond);
    println!("{:?}", sym_hash);
    println!("{:?}", subgraph_mult);

}

#![allow(non_snake_case)]

use itertools::izip;
use nauty_pet::graph::CanonGraph;
use petgraph::Undirected;
use std::collections::hash_map;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

fn _print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn cluster_to_iso_sym(
    cluster: &[usize],
    edge_list: &[Vec<usize>],
    iso_types: &[Vec<u8>],
    sym_types: &[Vec<u32>],
) -> (Vec<(usize, usize, u8)>, Vec<usize>) {
    // Takes cluster (sorted) and returns the corresponding ismorphic edge list and vertex type list
    let mut uncanon_list: Vec<(usize, usize, u8)> = vec![];
    let mut vertex_type: Vec<usize> = vec![];
    for (index, &vertex) in cluster.iter().enumerate() {
        vertex_type.push(0);
        for (&edge, &weight, &direction) in
            izip!(&edge_list[vertex], &iso_types[vertex], &sym_types[vertex])
        {
            let edge_position = cluster.iter().position(|&x| x == edge);
            match edge_position {
                Some(x) => {
                    vertex_type[index] += 2usize.pow(direction);
                    if x > index {
                        uncanon_list.push((index, x, weight));
                    }
                }
                _ => continue,
            }
        }
    }
    (uncanon_list, vertex_type)
}

fn cluster_to_iso(
    cluster: &[usize],
    edge_list: &[Vec<usize>],
    iso_types: &[Vec<u8>],
) -> Vec<(usize, usize, u8)> {
    // Takes cluster (unsorted) and returns the corresponding isomorphic edge list
    let mut uncanon_list: Vec<(usize, usize, u8)> = vec![];
    for (index, &vertex) in cluster.iter().enumerate() {
        for (&edge, &weight) in edge_list[vertex].iter().zip(iso_types[vertex].iter()) {
            let edge_position = cluster.iter().position(|&x| x == edge);
            match edge_position {
                Some(x) => {
                    if x > index {
                        uncanon_list.push((index, x, weight));
                    }
                }
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

fn add_cluster(
    mut cluster: Vec<usize>,
    lattice: &HashMap<usize, Vec<usize>>,
    edge_list: &[Vec<usize>],
    iso_types: &[Vec<u8>],
    sym_types: &[Vec<u32>],
    graph_multiplicity: &mut HashMap<usize, usize>,
    subgraph_multiplicity: &mut HashMap<usize, HashMap<usize, usize>>,
    graph_bond_info: &mut HashMap<usize, Vec<(usize, usize, u8)>>,
    sym_hash_set: &mut HashSet<usize>,
) {
    cluster.sort();

    let (cluster_iso_list, vertex_type_cluster) =
        cluster_to_iso_sym(&cluster, edge_list, iso_types, sym_types);

    let sym_hash = sym_to_hash(&vertex_type_cluster);

    if !sym_hash_set.contains(&sym_hash) {
        let iso_hash = iso_to_hash(&cluster_iso_list);

        match &mut graph_multiplicity.entry(iso_hash) {
            hash_map::Entry::Vacant(_) => {
                let mut lattice_clone = lattice.clone();
                lattice_clone.retain(|vertex, _| cluster.contains(vertex));

                sym_hash_set.insert(sym_hash);
                graph_multiplicity.insert(iso_hash, 1);
                graph_bond_info.insert(iso_hash, cluster_iso_list);

                let mut subgraph_func = |subcluster: Vec<usize>| {
                    add_subcluster(
                        subcluster,
                        edge_list,
                        iso_types,
                        iso_hash,
                        subgraph_multiplicity,
                    )
                };

                for size in 2..cluster.len() {
                    enumerate(&lattice_clone, size, &cluster, &mut subgraph_func);
                }
            }
            hash_map::Entry::Occupied(entry) => {
                sym_hash_set.insert(sym_hash);
                *entry.get_mut() += 1;
            }
        };
    };
}

fn add_subcluster(
    cluster: Vec<usize>,
    edge_list: &[Vec<usize>],
    iso_types: &[Vec<u8>],
    iso_hash: usize,
    subgraph_multiplicity: &mut HashMap<usize, HashMap<usize, usize>>,
) {
    let subcluster_iso_list = cluster_to_iso(&cluster, edge_list, iso_types);
    let sub_iso_hash = iso_to_hash(&subcluster_iso_list);
    subgraph_multiplicity
        .entry(iso_hash)
        .and_modify(|subgraph_info| {
            subgraph_info
                .entry(sub_iso_hash)
                .and_modify(|counter| *counter += 1)
                .or_insert(1);
        })
        .or_insert_with(|| HashMap::from([(sub_iso_hash, 1)]));
}

fn vsimple(
    edges: &HashMap<usize, Vec<usize>>,
    subgraph: &mut Vec<usize>,
    neighbors: &mut Vec<usize>,
    guarding_set: &HashSet<usize>,
    size: usize,
    graph_func: &mut dyn FnMut(Vec<usize>),
) -> bool {
    if subgraph.len() == size {
        graph_func(subgraph.clone());
        return true;
    };

    let mut has_int_leaf = false;

    let mut new_guarding_set = guarding_set.clone();

    while !neighbors.is_empty() {
        let neighbor = neighbors.pop().unwrap();

        if edges.contains_key(&neighbor) {
            subgraph.push(neighbor);

            let mut new_neighbors = neighbors.clone();

            for vertex in edges[&neighbor].iter() {
                if !subgraph.contains(vertex) & !new_guarding_set.contains(vertex) & !new_neighbors.contains(vertex) {
                    new_neighbors.push(*vertex);
                };
            };

            if vsimple(
                edges,
                subgraph,
                &mut new_neighbors,
                &new_guarding_set,
                size,
                graph_func,
            ) {
                subgraph.pop();
                has_int_leaf = true;
            } else {
                subgraph.pop();
                return has_int_leaf;
            };

            new_guarding_set.insert(neighbor);

            if (edges.len() - new_guarding_set.len()) < size {
                return has_int_leaf;
            };
        };
    }

    has_int_leaf
}

fn enumerate(
    edge_list: &HashMap<usize, Vec<usize>>,
    size: usize,
    starting_vertices: &Vec<usize>,
    mut graph_func: &mut dyn FnMut(Vec<usize>),
) {
    let mut guarding_set = HashSet::<usize>::new();

    for vertex in starting_vertices {
        let mut neighbors = edge_list[vertex]
            .iter()
            .cloned()
            .filter(|vertex| !guarding_set.contains(vertex))
            .collect();

        let mut starting_subgraph = vec![*vertex];
        vsimple(
            edge_list,
            &mut starting_subgraph,
            &mut neighbors,
            &guarding_set,
            size,
            &mut graph_func,
        );
        guarding_set.insert(*vertex);
    }
}

fn main() {
    use std::env;
    let args: Vec<_> = env::args().collect();

    let cluster_size: usize = args[1].parse().unwrap();
    let size: isize = cluster_size as isize;

    // Triangular Lattice
    //let directions = vec![(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)];
    //let weights = vec![1, 1, 1, 1, 1, 1];

    // Square Lattice
    let directions = vec![(1, 0), (0, 1), (-1, 0), (0, -1)];
    let weights = vec![1, 1, 1, 1];

    // Square Lattice nnn
    //let directions = vec![(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)];
    //let weights = vec![1, 1, 1, 1, 2, 2, 2, 2];

    // Triangular Lattice nnn
    //let directions = vec![(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1), (1, 2), (-1, 2), (1, -2), (2, 1), (-2, 1), (2, -1)];
    //let weights = vec![1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2];

    let conv = |(x, y): (isize, isize)| (x + (size * y)) as usize;

    let mut cluster_map = HashMap::<(isize, isize), Vec<((isize, isize), u8, u32)>>::new();

    for x in 0..size {
        for y in 0..size {
            let coord = (x, y);
            let info: Vec<((isize, isize), u8, u32)> = directions
                .clone()
                .into_iter()
                .enumerate()
                .map(|(index, (d1, d2))| ((d1 + x, d2 + y), weights[index], index as u32))
                .collect();
            cluster_map.insert(coord, info);
        }
    }
    let cluster_keys = cluster_map.clone();
    for (_, value) in cluster_map.iter_mut() {
        value.retain(|(coord, _, _)| cluster_keys.contains_key(coord));
    }

    let mut cluster = Vec::<usize>::new();
    let mut edge_list = Vec::<Vec<usize>>::new();
    let mut iso_types = Vec::<Vec<u8>>::new();
    let mut sym_types = Vec::<Vec<u32>>::new();
    cluster.resize(cluster_map.len(), 0);
    edge_list.resize(cluster_map.len(), Vec::new());
    iso_types.resize(cluster_map.len(), Vec::new());
    sym_types.resize(cluster_map.len(), Vec::new());

    for (key, value) in cluster_map.iter() {
        let converted = conv(*key);
        cluster[converted] = converted;
        for (edge, iso, sym) in value {
            edge_list[converted].push(conv(*edge));
            iso_types[converted].push(*iso);
            sym_types[converted].push(*sym);
        }
    }

    let start: usize = conv((size / 2, size / 2));
    println!("{}", start);

    let cluster_set: HashMap<usize, Vec<usize>> = cluster
        .iter()
        .cloned()
        .zip(edge_list.iter().cloned())
        .collect();

    let mut graph_mult = HashMap::<usize, usize>::new();
    let mut subgraph_mult = HashMap::<usize, HashMap<usize, usize>>::new();
    let mut graph_bond = HashMap::<usize, Vec<(usize, usize, u8)>>::new();
    let mut sym_hash = HashSet::<usize>::new();

    let mut graph_func = |cluster: Vec<usize>| {
        add_cluster(
            cluster,
            &cluster_set,
            &edge_list,
            &iso_types,
            &sym_types,
            &mut graph_mult,
            &mut subgraph_mult,
            &mut graph_bond,
            &mut sym_hash,
        )
    };

    use std::time::Instant;
    let now = Instant::now();

    enumerate(&cluster_set, cluster_size, &vec![start], &mut graph_func);

    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);

    let mut total = 0;
    let mut total_broken = 0;
    for (iso_hash, mult) in graph_mult.iter() {
        total += mult;
        total_broken += subgraph_mult[iso_hash].len();
    }
    println!("Isomorphically Distinct: {:?}", graph_mult.len());
    println!("Total: {:?}", total);
    println!("Total Broken Down: {:?}", total_broken);
    //println!("{:?}", graph_bond);
    //println!("{:?}", sym_hash);
    //println!("{:?}", subgraph_mult);
}

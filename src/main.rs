use itertools::izip;
use nauty_pet::graph::CanonGraph;
use petgraph::Undirected;
use std::collections::hash_map;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::fs::File;
use std::io::{Write};

fn cluster_to_uncanon_sym(
    cluster: &[usize],
    edge_list: &[Vec<usize>],
    iso_types: &[Vec<u8>],
    sym_types: &[Vec<u32>],
) -> (Vec<(usize, usize, u8)>, Vec<usize>) {
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

fn cluster_to_tr(
    cluster: &[usize]
) -> (Vec<isize>, usize) {
    // transformations in a closure: identity, 180 rotation, diagonal flip x = y, diagonal flip x = -y
    let transform = |x| (x, 136 - x, ((x * 17) - ((x / 17) * (288))), ((560 - (17 * x)) + (288 * (x / 17))));
    let cluster_orbit = cluster.iter().fold(vec![vec![], vec![], vec![], vec![]], |mut acc, x| {
        let temp = transform(*x as isize);
        acc[0].push(*x as isize);
        acc[1].push(temp.1);
        acc[2].push(temp.2);
        acc[3].push(temp.3);
        acc
    });
    let mut new_orbit: Vec<Vec<isize>> = vec![];
    for mut cl in cluster_orbit {
        cl.sort();
        let cl_min = cl[0];
        new_orbit.push(cl.iter().cloned().map(|x| x - cl_min).collect());
    }
    new_orbit.sort();

    let mut graph_hasher = hash_map::DefaultHasher::new();
    new_orbit.hash(&mut graph_hasher);
    (new_orbit[0].clone(), graph_hasher.finish() as usize)
}

//fn tr_to_hash(uncanon_list: &[(usize, usize, u8)]) -> usize {
//    let canon_graph = CanonGraph::<(), u8, Undirected, usize>::from_edges(uncanon_list);
//    let mut graph_hasher = hash_map::DefaultHasher::new();
//    canon_graph.hash(&mut graph_hasher);
//    graph_hasher.finish() as usize
//}

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
    subgraph_multiplicity: &mut HashMap<usize, HashMap<usize, (usize, usize)>>,
    graph_bond_info: &mut HashMap<usize, (Vec<(isize, isize)>, Vec<(usize, usize, u8)>)>,
    sym_hash_set: &mut HashSet<usize>,
) {
    cluster.sort();

    let (cluster_iso_list, vertex_type_cluster) =
        cluster_to_uncanon_sym(&cluster, edge_list, iso_types, sym_types);

    let sym_hash = sym_to_hash(&vertex_type_cluster);

    if !sym_hash_set.contains(&sym_hash) {
        let (form, iso_hash) = cluster_to_tr(&cluster);

        match &mut graph_multiplicity.entry(iso_hash) {
            hash_map::Entry::Vacant(_) => {
                let mut lattice_clone = lattice.clone();
                lattice_clone.retain(|vertex, _| cluster.contains(vertex));

                sym_hash_set.insert(sym_hash);
                graph_multiplicity.insert(iso_hash, 1);
                graph_bond_info.insert(iso_hash, (form.iter().map(|x| (x % 17, x / 17) ).collect(), cluster_iso_list ));

                let mut subgraph_func = |subcluster: Vec<usize>| {
                    add_subcluster(
                        subcluster,
                        edge_list,
                        iso_types,
                        iso_hash,
                        subgraph_multiplicity,
                    )
                };

                enumerate(&lattice_clone, 1, &cluster, &mut subgraph_func);
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
    subgraph_multiplicity: &mut HashMap<usize, HashMap<usize, (usize, usize)>>,
) {
    //let subcluster_iso_list = cluster_to_iso(&cluster, edge_list, iso_types);
    let (_, sub_iso_hash) = cluster_to_tr(&cluster);
    subgraph_multiplicity
        .entry(iso_hash)
        .and_modify(|subgraph_info| {
            subgraph_info
                .entry(sub_iso_hash)
                .and_modify(|(_order, counter)| *counter += 1)
                .or_insert((cluster.len(), 1));
        })
        .or_insert_with(|| HashMap::from([(sub_iso_hash, (cluster.len(), 1))]));
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
                if !subgraph.contains(vertex)
                    & !new_guarding_set.contains(vertex)
                    & !new_neighbors.contains(vertex)
                {
                    new_neighbors.push(*vertex);
                };
            }

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
    edges: &HashMap<usize, Vec<usize>>,
    size: usize,
    starting_vertices: &Vec<usize>,
    mut graph_func: &mut dyn FnMut(Vec<usize>),
) {
    let mut guarding_set = HashSet::<usize>::new();

    for vertex in starting_vertices {
        let mut neighbors = edges[vertex]
            .iter()
            .cloned()
            .filter(|neighbor| !guarding_set.contains(neighbor))
            .collect();

        let mut starting_subgraph = vec![*vertex];
        vsimple(
            edges,
            &mut starting_subgraph,
            &mut neighbors,
            &guarding_set,
            size,
            &mut graph_func,
        );
        guarding_set.insert(*vertex);
    }
}

fn gen_reg_lattice_2d(
    size: usize,
    directions: Vec<(isize, isize)>,
    weights: Vec<u8>,
) -> (usize, Vec<Vec<usize>>, Vec<Vec<u8>>, Vec<Vec<u32>>) {
    let mut cluster_map = HashMap::<(isize, isize), Vec<((isize, isize), u8, u32)>>::new();

    // Buffer size
    let size = 17;
    for x in 0..size as isize {
        for y in 0..size as isize {
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
    let cluster_keys: Vec<(isize, isize)> = cluster_map.keys().cloned().collect();

    for (_, value) in cluster_map.iter_mut() {
        value.retain(|(coord, _, _)| cluster_keys.contains(coord));
    }

    let mut edges = Vec::<Vec<usize>>::new();
    let mut iso_types = Vec::<Vec<u8>>::new();
    let mut sym_types = Vec::<Vec<u32>>::new();

    edges.resize(cluster_map.len(), Vec::new());
    iso_types.resize(cluster_map.len(), Vec::new());
    sym_types.resize(cluster_map.len(), Vec::new());

    // need to do again to get new keys or will end up with index out of bound errors
    let mut cluster_keys: Vec<(isize, isize)> = cluster_map.keys().cloned().collect();
    cluster_keys.sort();
    let conv = |key| cluster_keys.iter().position(|coord| coord == key).unwrap();

    for (key, value) in cluster_map.iter() {
        let converted = conv(key);
        for (edge, iso, sym) in value {
            edges[converted].push(conv(edge));
            iso_types[converted].push(*iso);
            sym_types[converted].push(*sym);
        }
    }

    let start: usize = conv(&(size as isize / 2, size as isize / 2));
    (start, edges, iso_types, sym_types)
}

fn main() -> std::io::Result<()>{
    use std::env;
    let args: Vec<_> = env::args().collect();

    let nlce_type: String = args[1].parse().unwrap();
    let cluster_size: usize = args[2].parse().unwrap();
    let nlce_directory = format!("./NLCE_Data/{}", nlce_type);

    let mut directions: Vec<(isize, isize)> = vec![];
    let mut weights: Vec<u8> = vec![];

    let options = ["triangle", "square", "square-next", "ani-triangle", "triangle-next"];
    match nlce_type.as_str() {
        "triangle" => {
    // Triangular Lattice
    directions = vec![(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)];
    weights = vec![1, 1, 1, 1, 1, 1];
        },
        "ani-triangle" => {
    // Square Off diagonal Lattice
    directions = vec![(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)];
    weights = vec![1, 2, 1, 1, 2, 1];
        },
        "square" => {
    // Square Lattice
    directions = vec![(1, 0), (0, 1), (-1, 0), (0, -1)];
    weights = vec![1, 1, 1, 1];
        },
        "square-next" => {
    // Square Lattice nnn
    directions = vec![(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)];
    weights = vec![1, 1, 1, 1, 2, 2, 2, 2];
        },
        "triangle-next" => {
    // Triangle Lattice nnn
    directions = vec![(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1), (-1, 1), (1, 2), (2, 1), (1, -1), (-2, -1), (-1, -2)];
    weights = vec![1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2];
        },
        _ => {
            println!("This is not a valid option, the current options are {:?}", options);
        }

    };

    let (start, edges, iso_types, sym_types) =
        gen_reg_lattice_2d(cluster_size, directions, weights);

    let cluster_vertices: Vec<usize> = (0..edges.len()).collect();
    let cluster_set: HashMap<usize, Vec<usize>> = cluster_vertices
        .iter()
        .cloned()
        .zip(edges.iter().cloned())
        .collect();

    let mut graph_mult = HashMap::<usize, usize>::new();
    let mut subgraph_mult = HashMap::<usize, HashMap<usize, (usize, usize)>>::new();
    let mut graph_bond = HashMap::<usize, (Vec<(isize, isize)>, Vec<(usize, usize, u8)>)>::new();
    let mut sym_hash = HashSet::<usize>::new();

    let mut graph_func = |cluster: Vec<usize>| {
        add_cluster(
            cluster,
            &cluster_set,
            &edges,
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


    let graph_mult_json = serde_json::to_string_pretty(&graph_mult)?;
    let subgraph_mult_json = serde_json::to_string_pretty(&subgraph_mult)?;
    let graph_bond_json = serde_json::to_string_pretty(&graph_bond)?;

    std::fs::create_dir_all(&nlce_directory).unwrap();

    let graph_mult_path = format!("{}/graph_mult_{}_{}.json", nlce_directory, nlce_type, cluster_size);
    let subgraph_mult_path = format!("{}/subgraph_mult_{}_{}.json", nlce_directory, nlce_type, cluster_size);
    let graph_bond_path = format!("{}/graph_bond_{}_{}.json", nlce_directory, nlce_type, cluster_size);

    let mut graph_mult_output = File::create(graph_mult_path)?;
    let mut subgraph_mult_output = File::create(subgraph_mult_path)?;
    let mut graph_bond_output = File::create(graph_bond_path)?;

    write!(graph_mult_output, "{}", graph_mult_json)?;
    write!(subgraph_mult_output, "{}", subgraph_mult_json)?;
    write!(graph_bond_output, "{}", graph_bond_json)?;

    Ok(())
}

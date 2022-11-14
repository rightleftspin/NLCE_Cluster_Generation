use petgraph::graphmap::UnGraphMap;

fn main() {
    let mut graph = UnGraphMap::<_, u8>::new();

    graph.add_node(0);
    graph.add_node(1);
    graph.add_edge(0, 1, 1);


    let mut graph1 = UnGraphMap::<_, u8>::new();

    graph1.add_node(0);
    graph1.add_node(1);
    graph1.add_edge(0, 1, 1);

}

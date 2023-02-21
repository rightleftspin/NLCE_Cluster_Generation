from sortedcontainers import SortedList
import sys

sys.setrecursionlimit(300)
# The main goal of this code is to provide a library for generating hashes from a specific graph
# in general, this
#

squareAdjecencyPossibilities = [(1, 0), (-1, 0), (0, 1), (0, -1)]
rotationLookupTable = [8, 4, 12, 1, 9, 5, 13, 2, 10, 6, 14, 3, 11, 7, 15]
flipLookupTable = [2, 1, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12, 14, 13, 15]

def rotation(vertexTypeGraph):
    rotatedVertexTypeGraph = ()
    for vertex in vertexTypeGraph:
        rotatedVertexTypeGraph += (rotationLookupTable[vertex - 1], )

    return(rotatedVertexTypeGraph)

def flip(vertexTypeGraph):
    flippedVertexTypeGraph = ()
    for vertex in vertexTypeGraph:
        flippedVertexTypeGraph += (flipLookupTable[vertex - 1], )

    return(flippedVertexTypeGraph)

transformationDict = { "r": [rotation], "f": [flip], "r2": [rotation, rotation], "r3": [rotation, rotation, rotation], "fr": [flip, rotation], "fr2": [flip, rotation, rotation], "fr3": [flip, rotation, rotation, rotation] }

def possibilityList(node, adjList):
    adjNodePossibilities = []
    for adjNode in adjList:
        adjNodePossibilities.append( tuple(coord + adj for coord, adj in zip(node, adjNode)) )

    return(adjNodePossibilities)

def isEnd(node, orderedGraph, nodeNumber):
    ind = orderedGraph.index(node)
    if ind < (len(orderedGraph) - 1):
        if node[0] < orderedGraph[ind + 1][0]:
            return((nodeNumber, 0))
        else:
            return((nodeNumber,))
    else:
        return((nodeNumber,))

def isometricNumbering(node, orderedGraph):
    nodeNumber = 0
    adjNodes = possibilityList(node, squareAdjecencyPossibilities)
    for adjNode in adjNodes:
        if adjNode in orderedGraph:
            nodeNumber += 1

    return(nodeNumber)

def symmetricNumbering(node, orderedGraph):
    nodeNumber = 0
    adjNodes = possibilityList(node, squareAdjecencyPossibilities)
    for index, adjNode in enumerate(adjNodes):
        if adjNode in orderedGraph:
            nodeNumber += 2 ** index

    return(nodeNumber)

def order(graph):
    # Takes in a graph in coordinate notation and orders it in a list with x values, then y values
    # this does the same as python's default sorted function
    return(sorted(graph))

def minimumVertexTypeNotation(graph, numberingFxn):
    # This function takes a graph and returns its minimum vertex-type tuple. The adjecencyOptions and numberingScheme
    # options detail the way in which the numberingFxn will convert a node into its corresponding
    # vertex type number.
    vertexTypeGraph = ()
    orderedGraph = order(graph)
    for node in orderedGraph:
        vertexTypeGraph += (numberingFxn(node, orderedGraph), )

    return(vertexTypeGraph)

def graphOrbitGenerator(vertexTypeGraph):
    # This needs to take a graph in vertex type form and return a list containing the orbit of the graph in vertex type form
    graphOrbit = [hash(tuple(sorted(vertexTypeGraph)))]
    orbitSize = 1
    for transformKey in transformationDict:
        newGraph = vertexTypeGraph
        for transform in transformationDict[transformKey]:
            newGraph = transform(newGraph)

        hashedTransformedGraph = hash(tuple(sorted(newGraph)))
        if hashedTransformedGraph not in graphOrbit:
            orbitSize += 1
            graphOrbit.append(hashedTransformedGraph)

    return(orbitSize, graphOrbit)


def checkGraphAlreadyExists(graph, isoNumberingFxn, symNumberingFxn, isoHashList, symHashList, graphDict):

    isoHash = hash(minimumVertexTypeNotation(graph, isoNumberingFxn))
    if isoHash not in isoHashList:
        symVertexTypeGraph = minimumVertexTypeNotation(graph, symNumberingFxn)
        symLen, symHashes = graphOrbitGenerator(vertexTypeGraph)

        graphDict[isoHash] = [graph, symLen]
        isoHashList.add(isoHash)
        symHashList.add(symHashes)
        return(isoHashList, symHashList, graphDict)
    else:
        symVertexTypeGraph = minimumVertexTypeNotation(graph, symNumberingFxn)
        if hash(symVertexTypeGraph) not in symHashList:
            symLen, symHashes = graphOrbitGenerator(vertexTypeGraph)

            graphDict[isoHash][1] += symLen
            symHashList.add(symHashes)

            return(isoHashList, symHashList, graphDict)

        else:
            return(isoHashList, symHashList, graphDict)



graphDictionary = {}
isomorphicHashList = SortedList()
symmetricHashList = SortedList()

testGraph1 = [(0,0), (1, 0), (0, 1), (1, 1), (2, 1), (1, 2)]
testGraph2 = []
gCounter = set()

def symHash(graph):
    minGraph = min(graph)
    shiftedGraph = frozenset(map(lambda x: (x[0] - minGraph[0], x[1] - minGraph[1]), graph))
    return(hash(shiftedGraph))

def enumerateGraph(graph, size):
    guardingSet = set()
    #for vertex in graph:
    vertex = (5, 5)
    neighbors = graph[vertex]
    vSimple(graph, {vertex}, neighbors.difference(guardingSet), guardingSet, size)
    guardingSet = guardingSet | {vertex}

def vSimple(graph, subgraph, neighbors, guardingSet, size):
    if len(subgraph) == size:
        global gCounter
        gCounter.add(hash(minimumVertexTypeNotation(subgraph, symmetricNumbering)))
        return(True)
    
    hasIntLeaf = False
    for neighbor in neighbors:
        newSubgraph = subgraph | {neighbor}
        addNeighbors = graph[neighbor]
        newNeighbors = neighbors.difference({neighbor}) | (addNeighbors.difference(newSubgraph).difference(guardingSet))

        if vSimple(graph, newSubgraph, newNeighbors, guardingSet, size):
            hasIntLeaf = True
        else:
            break

        guardingSet = guardingSet | {neighbor}
        if (len(graph) - len(guardingSet)) < size:
            break

    return(hasIntLeaf)


#Generate test 11 by 11 graph
size = 11
testVert = set()
testGraph = {}
for x in range(size):
    for y in range(size):
        testVert.add((x, y))
        testGraph[(x, y)] = set()


squareAdjecencyPossibilities = [(1, 0), (-1, 0), (0, 1), (0, -1)]
for vert in testVert:
    for adj in squareAdjecencyPossibilities:
        edge = (vert[0] + adj[0], vert[1] + adj[1])
        if edge in testVert:
            testGraph[vert].add(edge)

#print(testGraph)
    
#testGraph = {(0, 0): {(1, 0), (0, 1)}, (1, 0): {(0, 0), (1, 1)}, (1, 1):{(1, 0), (0, 1)}, (0, 1):{(1, 1), (0, 0)}}
import time

start = time.time()
enumerateGraph(testGraph, 7)
print(time.time() - start)
print(len(gCounter))

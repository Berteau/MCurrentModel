import networkx as nx
from pylab import *

class NetPlotter:
    def __init__(self, connections = None):
        self.G = nx.DiGraph()
        if connections == None:
            self.connections = []
        else:
            self.connections = connections

    def plot(self):
        weights = [connection.weight for connection in self.connections]
        maxWeight = max(weights)
        tempEdges = [(connection.source.name, connection.target.name, connection.weight/maxWeight) for connection in self.connections]
        self.G.add_weighted_edges_from(tempEdges)

        pos = nx.spring_layout(self.G)
        nx.draw_networkx_nodes(self.G, pos, node_size=100)
        nx.draw_networkx_labels(self.G, pos)

        excite = [(u, v) for (u, v, d) in self.G.edges(data=True) if d['weight'] > 0]
        inhibit = [(u, v) for (u, v, d) in self.G.edges(data=True) if d['weight'] < 0]

        nx.draw_networkx_edges(self.G, pos, edgelist=excite, edge_color='g', arrows=True)
        nx.draw_networkx_edges(self.G, pos, edgelist=inhibit, edge_color='r', arrows=True)

        show()





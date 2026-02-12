"""
Patch to fix GraphML export issue in tissue_simulator.

The issue: NetworkX's write_graphml() doesn't support numpy arrays as node attributes.
The solution: Convert numpy arrays to lists before exporting.

Apply this fix to tissue_simulator/spatial_analysis.py
"""

# OLD CODE (around line 452-468):
OLD_CODE = '''
    def export_network(self, filename: str, format: str = "graphml"):
        """
        Export the network to a file.
        
        Args:
            filename: Output filename
            format: "graphml", "gexf", "gml", or "edgelist"
        """
        if self.graph is None:
            raise ValueError("Network not built.")
        
        if format == "graphml":
            nx.write_graphml(self.graph, filename)
        elif format == "gexf":
            nx.write_gexf(self.graph, filename)
        elif format == "gml":
            nx.write_gml(self.graph, filename)
        elif format == "edgelist":
            nx.write_edgelist(self.graph, filename, data=['weight', 'distance'])
        else:
            raise ValueError(f"Unknown format: {format}")
'''

# NEW CODE (REPLACEMENT):
NEW_CODE = '''
    def export_network(self, filename: str, format: str = "graphml"):
        """
        Export the network to a file.
        
        Args:
            filename: Output filename
            format: "graphml", "gexf", "gml", or "edgelist"
        """
        if self.graph is None:
            raise ValueError("Network not built.")
        
        if format == "graphml":
            # GraphML doesn't support numpy arrays - convert to lists
            graph_copy = self.graph.copy()
            for node, data in graph_copy.nodes(data=True):
                for key, value in data.items():
                    if isinstance(value, np.ndarray):
                        data[key] = value.tolist()
            nx.write_graphml(graph_copy, filename)
        elif format == "gexf":
            nx.write_gexf(self.graph, filename)
        elif format == "gml":
            nx.write_gml(self.graph, filename)
        elif format == "edgelist":
            nx.write_edgelist(self.graph, filename, data=['weight', 'distance'])
        else:
            raise ValueError(f"Unknown format: {format}")
'''

print("=" * 80)
print("FIX FOR GRAPHML EXPORT ERROR")
print("=" * 80)
print("\nError: GraphML writer does not support <class 'numpy.ndarray'> as data values")
print("\nSOLUTION:")
print("-" * 80)
print(NEW_CODE)
print("-" * 80)
print("\nINSTRUCTIONS:")
print("1. Open: tissue_simulator/tissue_simulator/spatial_analysis.py")
print("2. Find the export_network() method (around line 452)")
print("3. Replace the graphml section with the code above")
print("\nThe fix creates a copy of the graph and converts all numpy arrays")
print("to Python lists before exporting to GraphML format.")
print("=" * 80)
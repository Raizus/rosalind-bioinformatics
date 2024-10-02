from typing import Any
import graphviz
import networkx as nx

def draw_graph(dot: graphviz.Graph, graph: nx.Graph):
    # Add nodes for molecules and components
    for node_id, node_data in graph.nodes(data=True):
        label_data = node_data['label']
        molecule_name = label_data[0]
        component_name = label_data[1]
        state = label_data[2]
        bond = label_data[3]
        label = ""
        shape = ""
        style = ""

        # Molecule node
        if component_name is None:
            label = f"({molecule_name}, {node_id[0]})"
            shape = 'rectangle'
        else:
            # Component node
            molecule_idx, component_idx = node_id
            label = f"({component_name}, {molecule_idx}, {component_idx})"
            if state:
                label += f"~{state}"
            shape = 'ellipse'

            # Determine border style
            if bond == '?':
                style = 'dashed'
            elif bond:
                style = 'dotted'
            else:
                style = 'solid'

        dot.node(str(node_id), label=label, shape=shape, style=style)

    # Add edges between molecules and components
    for edge_start, edge_end, edge_data in graph.edges(data=True):
        bond = edge_data.get('bond', None)
        style = 'solid'
        label = ""

        # If it's a bond between two components
        if bond:
            style = 'dashed'
            label = bond

        dot.edge(str(edge_start), str(edge_end), style=style, label=label)

    return dot


def format_data_into_lines(data: list[tuple[Any,...]]):
    lines: list[str] = []

    col_widths = [max(len(str(row[i])) for row in data)
                    for i in range(len(data[0]))]
    for row in data:
        line = "\t".join(f"{str(item):<{col_widths[i]}}" 
                        for i, item in enumerate(row))
        lines.append(f"\t{line}")

    return lines


def lowest_missing_positive_integer(string_list: list[str]) -> str:
    # Convert valid strings to integers
    integers = set(int(x) for x in string_list if x.isdigit())

    # Start from 1 and find the first missing positive integer
    i = 1
    while i in integers:
        i += 1

    return str(i)

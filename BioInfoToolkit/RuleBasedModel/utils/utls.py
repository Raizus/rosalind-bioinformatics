
import ast
import operator as op
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


# Supported operators
ALLOWED_OPERATORS = {
    ast.Add: op.add,
    ast.Sub: op.sub,
    ast.Mult: op.mul,
    ast.Div: op.truediv,
    ast.Pow: op.pow,
    ast.BitXor: op.xor,
    ast.USub: op.neg
}


def eval_expr(expr: str, variables: dict[str, Any]):
    """
    Safely evaluate an algebraic expression with variables.
    
    Args:
    - expr (str): The expression to evaluate.
    - variables (dict): A dictionary of variable names and their values.
    
    Returns:
    - result: The evaluated result of the expression.
    """

    def eval_node(node):
        if isinstance(node, ast.Constant):  # <number>
            return node.n

        if isinstance(node, ast.BinOp):  # <left> <operator> <right>
            left = eval_node(node.left)
            right = eval_node(node.right)

            op_type = type(node.op)
            if op_type in ALLOWED_OPERATORS:
                return ALLOWED_OPERATORS[op_type](left, right)

            raise TypeError(f"Unsupported binary operator: {op_type}")

        if isinstance(node, ast.UnaryOp):  # - <operand> e.g., -1
            operand = eval_node(node.operand)
            op_type = type(node.op)
            if op_type in ALLOWED_OPERATORS:
                return ALLOWED_OPERATORS[op_type](operand)

            raise TypeError(f"Unsupported unary operator: {op_type}")

        if isinstance(node, ast.Name):  # variable
            if node.id in variables:
                return variables[node.id]
            raise NameError(f"Use of undefined variable: {node.id}")

        raise TypeError(f"Unsupported type: {type(node)}")

    # Parse expression
    parsed_expr = ast.parse(expr, mode='eval').body

    # Evaluate parsed AST
    return eval_node(parsed_expr)

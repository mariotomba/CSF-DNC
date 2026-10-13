"""FINITE_EXACT regression checks for equal-order leaf components.

Run from the repository root with:
    sage -python -m unittest discover -s tests -v

Only notebook definitions are executed; the research example is not run.
"""

import contextlib
import io
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from sage.all import Graph, Integer, Partitions, QQ, SymmetricFunctions, binomial, graphs
from sage.plot.graphics import Graphics
from sage.repl.preparse import preparse


NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks" / "DNCfundamentals.ipynb"


def load_definitions():
    namespace = {"__name__": "dnc_notebook_definitions"}
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    for index, cell in enumerate(notebook["cells"]):
        source = "".join(cell["source"])
        if cell["cell_type"] == "markdown" and source.lstrip().startswith("## Example"):
            break
        if cell["cell_type"] == "code":
            exec(compile(preparse(source), f"{NOTEBOOK}:cell {index}", "exec"), namespace)
    return namespace


DNC = load_definitions()


def canonical_key(tree):
    return tree.canonical_label().graph6_string()


def leaf_component_orders(tree):
    """Independent definition: remove edges with two original non-leaf ends."""
    original_degrees = dict(tree.degree_iterator(labels=True))
    forest = tree.copy()
    forest.delete_edges(
        (u, v)
        for u, v in tree.edges(labels=False)
        if original_degrees[u] > 1 and original_degrees[v] > 1
    )
    return sorted(len(component) for component in forest.connected_components())


def power_sum_csf(tree, p):
    """Stanley's spanning-subgraph formula, independently of the DNC recurrence."""
    vertices = tree.vertices(sort=False)
    edges = tree.edges(labels=False)
    result = p.zero()
    for mask in range(1 << len(edges)):
        forest = Graph()
        forest.add_vertices(vertices)
        selected = [edge for index, edge in enumerate(edges) if mask & (1 << index)]
        forest.add_edges(selected)
        partition = sorted((len(c) for c in forest.connected_components()), reverse=True)
        result += (-1) ** len(selected) * p[partition]
    return result


def star_vector_in_power_sums(vector, n, p):
    # The connected-tree output uses st[1] = p[1]. The helper's signed singleton
    # base cases implement the signs of its all-plus recurrence, not a new basis.
    stars = {1: p[1]}
    for size in range(2, n + 1):
        stars[size] = sum(
            (
                (-1) ** kept * binomial(size - 1, kept)
                * p[[kept + 1] + [1] * (size - kept - 1)]
                for kept in range(size)
            ),
            p.zero(),
        )
    result = p.zero()
    for coefficient, partition in zip(vector, Partitions(n)):
        if coefficient:
            term = p.one()
            for part in partition:
                term *= stars[part]
            result += Integer(coefficient) * term
    return result


def capture_display(function, *args, **kwargs):
    output = io.StringIO()
    before_plots = []
    with contextlib.redirect_stdout(output), patch.object(
        Graphics, "show", side_effect=lambda *a, **k: before_plots.append(output.getvalue())
    ) as show:
        result = function(*args, **kwargs)
    return result, output.getvalue(), show.call_count, before_plots


class UniformLeafComponentTests(unittest.TestCase):
    def test_matches_exhaustive_unlabeled_enumeration_through_ten_vertices(self):
        for n in range(2, 11):
            expected_by_order = {l: set() for l in range(2, n + 1)}
            for tree in graphs.trees(n):
                orders = leaf_component_orders(tree)
                if orders and orders[0] >= 2 and len(set(orders)) == 1:
                    expected_by_order[orders[0]].add(canonical_key(tree))
            for l, expected in expected_by_order.items():
                with self.subTest(n=n, l=l):
                    produced = list(DNC["trees_uniform_leaf_components"](n, l))
                    actual = {canonical_key(tree) for tree in produced}
                    self.assertEqual(actual, expected)
                    self.assertEqual(len(produced), len(actual), "Duplicate isomorphism class")
                    for tree in produced:
                        self.assertTrue(tree.is_tree())
                        self.assertEqual(tree.order(), n)
                        self.assertEqual(leaf_component_orders(tree), [l] * (n // l))

    def test_single_leaf_component_is_the_star_including_k2(self):
        for n in range(2, 9):
            with self.subTest(n=n):
                trees = list(DNC["trees_uniform_leaf_components"](n, n))
                self.assertEqual(len(trees), 1)
                self.assertTrue(trees[0].is_isomorphic(graphs.StarGraph(n - 1)))
                vector = DNC["CSF_tree"](trees[0])
                self.assertEqual(list(vector), [1] + [0] * (len(Partitions(n)) - 1))

    def test_accepts_python_sage_and_numpy_integers(self):
        for constructor in (int, Integer, np.int64):
            with self.subTest(integer_type=constructor):
                trees = list(DNC["trees_uniform_leaf_components"](constructor(4), constructor(2)))
                self.assertEqual(len(trees), 1)
                self.assertTrue(trees[0].is_isomorphic(graphs.PathGraph(4)))

    def test_invalid_inputs_are_rejected_by_both_interfaces(self):
        invalid = (-3, 0, 1, 2.0, 2.5, "2", None, True, False, np.bool_(True))
        for value in invalid:
            for args in ((value, 2), (4, value)):
                with self.subTest(args=args):
                    with self.assertRaises(ValueError):
                        list(DNC["trees_uniform_leaf_components"](*args))
                    with self.assertRaises(ValueError):
                        capture_display(DNC["CSF_uniform_leaf_components"], *args)

    def test_no_solutions_have_no_plots_and_a_clear_message(self):
        for n, l in ((3, 2), (7, 3), (2, 3)):
            with self.subTest(n=n, l=l):
                self.assertEqual(list(DNC["trees_uniform_leaf_components"](n, l)), [])
                result, output, plots, _ = capture_display(DNC["CSF_uniform_leaf_components"], n, l)
                self.assertIsNone(result)
                self.assertEqual(plots, 0)
                self.assertIn("no trees", output.lower())
                self.assertIn(str(n), output)
                self.assertIn(str(l), output)

    def test_display_matches_existing_csf_tree_output_and_plot_order(self):
        for n, l in ((2, 2), (5, 5), (8, 2), (9, 3)):
            with self.subTest(n=n, l=l):
                trees = list(DNC["trees_uniform_leaf_components"](n, l))

                def existing_display():
                    for tree in trees:
                        DNC["CSF_tree"](tree, prnt=1)

                _, expected, expected_plots, expected_before = capture_display(existing_display)
                result, actual, actual_plots, actual_before = capture_display(
                    DNC["CSF_uniform_leaf_components"], n, l
                )
                self.assertIsNone(result)
                self.assertEqual(actual, expected)
                self.assertEqual(actual_plots, expected_plots)
                self.assertEqual(actual_plots, len(trees))
                self.assertEqual(actual_before, expected_before)

    def test_star_expansions_equal_independent_power_sum_csf(self):
        p = SymmetricFunctions(QQ).p()
        for n, l in ((2, 2), (4, 2), (5, 5), (6, 2), (8, 2), (9, 3)):
            for tree in DNC["trees_uniform_leaf_components"](n, l):
                with self.subTest(n=n, l=l, tree=canonical_key(tree)):
                    vector = DNC["CSF_tree"](tree)
                    actual = star_vector_in_power_sums(vector, n, p)
                    self.assertEqual(actual, power_sum_csf(tree, p))


if __name__ == "__main__":
    unittest.main()

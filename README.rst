cgnav module
============

The `cgnav` module builds a *pattern concept lattice* [GK01]_ (with *abstract concept graphs* [Wil97]_ as intents) from a *power context family* [Wil97]_.

Sample input files ``resources/family_relations.pcf`` and ``resources/circles.pcf`` are provided.
They contain power context families, from which the concept lattices presented
in [Koe13]_ (Fig. 2, :math:`\mathfrak{C}_{\mathcal{F}}[\{x\}]` only) and [Koe16a]_ (Fig. 7) are built.

**Usage:** The following code example generates the concept lattice in the example below.

.. code-block:: python

  >>> from cgnav import *
  >>> pcf = PCFFile("family_relations.pcf").parse()
  >>> lattice = PatternStructure(IntensionGraph,pcf).build()
  >>> top_concept = lattice[-1]

Concept intents, more generally intension graphs, also power context families, and a bit more, can be displayed in Jupyter Notebook (but
`pygraphviz <https://pygraphviz.github.io/>`_ must be installed).
An example notebook ``cstiw16.ipynb``, originally presented at the *Conceptual Structures Tools & Interoperability Workshop* (CSTIW 2016) in Annecy, France,
and the associated publication ``resources/cstiw16.pdf`` [Koe16b]_, are provided with the module.

**Example:**
The power context family ``family_relations.pcf`` is represented by the two cross tables below on the left. It describes five persons: A(nne), B(ob), C(hris),
D(ora) and E(mily). The second table states that Anne is the mother of Bob and Chris, and Bob is the father of Dora and Emily.
The concept lattice :math:`\mathfrak{C}_{\mathcal{F}}[\{x\}]`, reproduced below on the right, is a hierarchy of nine concepts; these are abstractions of
the data in the power context family. Each concept is described by a pattern (a labeled graph with a yellow node) and a set of persons (the examples of the
concept); the downward lines connect each concept with its subconcepts. The pattern is called the *concept intent*,
and the set is called the *concept extent*.

.. list-table::
   :widths: 40 60

   * - .. raw:: html

          <style>
          .pcf tr:nth-child(even) {background-color: #f5f5f5;}
          .pcf tr td {text-align: right;}
          </style>

          <table class="pcf" style='display:inline-block'>
          <tr><td></td><td>male</td><td>female</td></tr>
          <tr><td>A</td><td></td><td>x</td></tr>
          <tr><td>B</td><td>x</td><td></td></tr>
          <tr><td>C</td><td>x</td><td></td></tr>
          <tr><td>D</td><td></td><td>x</td></tr>
          <tr><td>E</td><td></td><td>x</td></tr>
          </table>&nbsp;&nbsp;<table class="pcf" style='display:inline-block'>
          <tr><td></td></tr>
          </table>&nbsp;&nbsp;<table class="pcf" style='display:inline-block'>
          <tr><td></td><td>father</td><td>mother</td><td>parent</td></tr>
          <tr><td>(A, B)</td><td></td><td>x</td><td>x</td></tr>
          <tr><td>(A, C)</td><td></td><td>x</td><td>x</td></tr>
          <tr><td>(B, D)</td><td>x</td><td></td><td>x</td></tr>
          <tr><td>(B, E)</td><td>x</td><td></td><td>x</td></tr>
          </table>

     - .. image:: resources/cube.svg

The most general concept, at the very top, describes all *persons*; it's pattern consists of a single node
(in this example, pattern nodes represent persons). Its subconcept on the right describes all *females*,
as indicated by the pattern node's label. The top concept is connected to two further subconcepts; their concept intents differ
only in the position of the yellow node. The p-labeled pattern edge represents a parent relation, pointing from the parent to the child;
the yellow node designates the subject of the description; so the left subconcept describes all *children*, and the middle subconcept describes
all *parents*. Naturally, the concept just below *children* and *females* should represent all *daughters*;
and indeed, its extent contains exactly the daughters, Dora and Emily; but the concept intent has further information,
it reveals that Dora and Emily are in fact *granddaughters*. The concept intents describe everything (representable in a pattern)
that the individuals in the concept extent have in common.

**Background:** 
A classical paper, *Lattice Model of Browsable Data Spaces* [GSG86]_, uses concept lattices for a keyword-based search of documents.
The concepts are places in an information space; at each place, you can find a collection of objects (the concept extent),
and a description of these objects (the concept intent). The lines between the concepts (as shown in the diagram) are paths
that take you from one place to another. *Concept Lattices of a Relational Structure* [Koe13]_ mathematically defines such information spaces e.g. for
relational databases (the above concept lattice is an example from that paper). *A Database Browser based on Pattern Concepts* [KS13]_
(`slides <https://www.hse.ru/data/2013/05/23/1299138754/Kotters%20Shmidt%20A%20Database%20Browser%20based%20on%20Pattern.pdf>`_)
offers a command-line interface to explore the concept lattice, but presenting concept intents as formulas on the command-line was not
practicable. The `cgnav` module is a re-implementation of the core database browser that

* supports visual display e.g. of concept intents in the Jupyter Notebook,
* is based on a re-formulation [Koe16a]_ of *Concept Lattices of a Relational Structure* in the notational framework of an inspirational paper
  by Rudolf Wille [Wil97]_.

**Limitations:**
The implementation of the lattice building algorithm can handle only small power context families, like the two provided as examples.
Only a certain kind of concept lattices is computed, where pattern intents describe single objects (i.e. have a single yellow node).

The idea of concept lattices as navigation spaces is appealingly simple, but especially in the given setting,
a more refined approach is required:

* concept intents are in general too complex to be understood
* concept intents do not reflect the user's intention in an exploratory search
* limiting movement to the edges excludes meaningful navigation options

The `dbnav <https://github.com/koetters/dbnav>`_ project addresses these problems,
but supports no lattice building, and uses a minimalistic version of concept intents.

**References:**

.. [GK01] Bernhard Ganter and Sergei Kuznetsov, *Pattern Structures and Their Projections*, Proceedings of ICCS 2001, LNCS 2120, 129-142 (2001)
.. [GSG86] Robert Godin, Eugene Saunders and Jan Gecsei, *Lattice Model of Browsable Data Spaces*, Inf. Sci. 40(2), 89-116 (1986)
.. [Koe13] Jens Kötters, *Concept Lattices of a Relational Structure*, Proceedings of ICCS 2013, LNCS 7735, 301-310 (2013)
.. [KS13] Jens Kötters and Heinz W. Schmidt, *A Database Browser based on Pattern Concepts*, Proceedings of FCAIR 2013, CEUR 977, 47-56 (2013)
.. [Koe16a] Jens Kötters, *Intension Graphs as Patterns over Power Context Families*, Proceedings of CLA 2016, CEUR 1624, 203-216 (2016)
.. [Koe16b] Jens Kötters, *A Python Library for FCA with Conjunctive Queries*, Proceedings of CSTIW 2016, CEUR 1637, 55-62 (2016)
.. [Wil97] Rudolf Wille, *Conceptual Graphs and Formal Concept Analysis*, Proceedings of ICCS 1997, LNCS 1257, 290-303 (1997)


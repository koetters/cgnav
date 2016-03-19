CGNav
=====
Constructs a pattern concept lattice (with abstract concept graphs as intents) from a power context family.
The file "familyrelations.pcf" can be used as a sample input file. The file format is an extension of
Burmeister's cxt format (although there may be small differences as to what parses correctly).

Running "python cgnav.py <input.pcf>" builds the lattice and starts an interactive lattice navigation session
in the lattice's top concept. Use keys up(go up),down(go down),u(show upper cover),l(show lower cover).

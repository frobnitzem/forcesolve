# A standard, vanilla forcefield using terms from the 1980s.
terms = [
    Term("pbond",     PolyBond,    Conn(1,2) ),
    Term("pangle",    PolyAngle,   Conn(1,2,3) ),
    Term("pub",       PolyBond,    Conn(1,3) ),
    Term("ptor",      PolyTorsion, Conn(1,2,3,4) ),
    Term("pimproper", PolyImprop,  OOP() ),
    Term("ljpair_1,4", LJPair, Conn(1,4) ),
    PairTerm("ljpair_5+",  LJPair, Conn(1,2) | Conn(1,3) | Conn(1,4) )
]


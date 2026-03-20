"""Physical constants and shared lookup tables used throughout the LaREST pipeline.

Unit conversion factors convert values produced by external QM programs (xTB,
CENSO) into SI-consistent J/mol units before they are stored or compared.

SMARTS patterns in ``MONOMER_GROUPS`` and ``INITIATOR_GROUPS`` define the
ring-opening functional groups that LaREST recognises when building polymer
chains.
"""

HARTTREE_TO_JMOL = 2625499.63948
KCALMOL_TO_JMOL = 4184.0
CALMOL_TO_JMOL = 4.184
MONOMER_GROUPS: dict[str, str] = {
    "CC": "[O;R]-[C;R](=[O;!R])-[O;R]",
    "CtC": "[O;R]-[C;R](=[O;!R])-[S;R]",
    "CdtC": "[S;R]-[C;R](=[O;!R])-[S;R]",
    "CtnC": "[O;R]-[C;R](=[S;!R])-[O;R]",
    "CX": "[O;R]-[C;R](=[S;!R])-[S;R]",
    "CtX": "[S;R]-[C;R](=[S;!R])-[S;R]",
    "L": "[C,c;R]-[C;R](=[O;!R])-[O;R]",
    "tL": "[C,c;R]-[C;R](=[O;!R])-[S;R]",
    "tnL": "[C,c;R]-[C;R](=[S;!R])-[O;R]",
    "dtL": "[C,c;R]-[C;R](=[S;!R])-[S;R]",
    "oA": "[O;R]-[C;R](=[O;!R])-[N;R]",
    "Lm": "[C,c;R]-[C;R](=[O;!R])-[N;R]",
}
# TODO: include additional initiator groups for further customisation
INITIATOR_GROUPS: dict[str, str] = {"R-OH": "[O;H][H]"}
CREST_ENTROPY_OUTPUT_PARAMS: list[str] = [
    "S_conf",
    "S_rrho",
    "S_total",
]
THERMODYNAMIC_PARAMS: list[str] = [
    "H",  # E taken as H for CENSO
    "S",
    "G",
]

CENSO_SECTIONS: list[str] = [
    "censo_prescreening",
    "censo_screening",
    "censo_optimization",
    "censo_refinement",
]
PIPELINE_SECTIONS: list[str] = [
    "rdkit",
    "crest",
    "censo_prescreening",
    "censo_screening",
    "censo_optimization",
    "censo_refinement",
]
ENTHALPY_PLOTTING_SECTIONS: list[str] = [
    "rdkit",
    "crest",
    "censo_refinement",
]
ENTROPY_PLOTTING_SECTIONS: list[str] = [
    "rdkit",
    "crest",
    "censo_refinement",
    "censo_corrected",
]

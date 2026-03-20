"""LaREST: Lactone Ring-opening Energetics Sorting Tool.

A computational chemistry pipeline for calculating thermodynamic parameters
(H, S, G) of lactone ring-opening polymerization reactions. The pipeline
chains together RDKit, xTB, CREST, CENSO, and ORCA to generate and refine
conformer ensembles, ultimately predicting ring-opening energetics.

Authors
-------
Ryan Reese, Jacky Zhang, Hannah Wilding, Alexander Ganose,
Victor Posligua-Hernandez, Charles Romain

Department of Chemistry, Imperial College London
"""

LAREST_HEADER: str = """
       ===========================================================
       |                                                         |
       |                        LaREST                           |
       |                                                         |
       |       Lactone Ring-opening Energetics Sorting Tool      |
       |                                                         |
       |         Ryan Reese, Jacky Zhang, Hannah Wilding,        |
       |        Alexander Ganose, Victor Posligua-Hernandez      |
       |                    Charles Romain                       |
       |                                                         |
       |    Department of Chemistry, Imperial College London     |
       |                                                         |
       ===========================================================
"""

#!/usr/bin/env python3

# ················································································· #

# Imports
import argparse
import itertools
import pandas as pd
import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley


# ················································································· #

def normalize_dict(dict: dict,
                   final_min: int,
                   final_max: int):
    """
    Utility function to normalize the values of a dictionary to meet a specified minimum and maximum value.
    
    :param dict: The dictionary of numeric values to normalize.
    :param final_min: The minimum value of the dictionary after normalization.
    :param final_max: The maximum value of the dictionary after normalization.
    :return: A dictionary normalized with the specified minimum/maximum values.
    """
    # Determining initial min/max
    init_min = min(dict.values())
    init_max = max(dict.values())
    
    # Normalizing to min=0, max=1
    norm_dict = {key: (value - init_min) / (init_max - init_min) for key, value in dict.items()}
    
    # Scaling to desired min/max values
    return {key: (value * (final_max - final_min)) + final_min for key, value in norm_dict.items()}


# ················································································· #

# Hydrophobicity scales (Kyte & Doolittle #TODO: Verify)
HYDROPHOBICITY = {
    "ALA": 1.8,  "ARG": -4.5, "ASN": -3.5, "ASP": -3.5,
    "CYS": 2.5,  "GLN": -3.5, "GLU": -3.5, "GLY": -0.4,
    "HIS": -3.2, "ILE": 4.5,  "LEU": 3.8,  "LYS": -3.9,
    "MET": 1.9,  "PHE": 2.8,  "PRO": -1.6, "SER": -0.8,
    "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2
}

NORMALIZED_HYDROPHOBICITY = normalize_dict(HYDROPHOBICITY, final_min=1, final_max=2)

# Charge definitions
CHARGE = {'ASP': -1, 'GLU': -1, 'LYS': 1, 'ARG': 1, 'HIS': 0.1}

# IMGT CDR definitions (Residues have IMGT numbering in ABodyBuilder2-generated PDB)
IMGT_LC_CDRS = [*range(27, 38+1), *range(56, 65+1), *range(105, 117+1)]
IMGT_HC_CDRS = [*range(27, 38+1), *range(56, 65+1), *range(105, 117+1)]

# Defining IMGT anchor residue positions as residues immediately neighbouring IMGT defined CDR positions
IMGT_LC_ANCHORS = (set(pos + 1 for pos in IMGT_LC_CDRS) | set(pos - 1 for pos in IMGT_LC_CDRS)) - set(IMGT_LC_CDRS)
IMGT_HC_ANCHORS = (set(pos + 1 for pos in IMGT_HC_CDRS) | set(pos - 1 for pos in IMGT_HC_CDRS)) - set(IMGT_HC_CDRS)

# SASA of Ala-Arg-Ala tripeptide with the Shrake and Rupley algorithm
ARA_SASA = 241 #TODO: Verify

# Backbone atom names in PDB
backbone_atoms = {"N", "CA", "C", "O"}


# ················································································· #

def load_pdb(pdb_path: str):
    """
    Takes the path to a PDB structure,
    returns the first model in the structure.
    
    :param pdb_path: The path to a PDB structure to parse with ``Bio.PDB.PDBParser``.
    :return model: A PDB.Model.Model object of the PDB file.
    """
    # Parsing structure and returning first model
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.replace('.pdb',''), pdb_path)
    model = structure[0]
    return model


# ················································································· #

def cdr(model: PDB.Model.Model) -> set[PDB.Residue.Residue]:
    """
    Takes a PDB model with residue IDs corresponding to IMGT numbering, returns the CDR residues of the structure.
    A PDB of a model from ABodyBuilder2 is indexed according to IMGT numbering.
    
    :param model: A IMGT-numbered PDB model.
    :return cdr_residues: A set of the CDR residues of the model.
    """
    # Initializing set
    cdr_residues = set()
    
    # Determining CDR residues in light and heavy chain by whether their index matches an IMGT-defined CDR-position
    for residue in model['L']:
        if residue.id[1] in IMGT_LC_CDRS:
            cdr_residues.add(residue)
    for residue in model['H']:
        if residue.id[1] in IMGT_HC_CDRS:
            cdr_residues.add(residue)
            
    return cdr_residues


def anchor(model: PDB.Model.Model) -> set[PDB.Residue.Residue]:
    """
    Takes a PDB model with residue IDs corresponding to IMGT numbering, returns the anchor residues of the structure
    (i.e., residues immediately flanking the CDRs).
    A PDB of a model from ABodyBuilder2 is indexed according to IMGT numbering.
    
    :param model: A IMGT-numbered PDB model.
    :return anchor_residues: A set of the anchor residues of the model.
    """
    # Initializing set
    anchor_residues = set()
    
    # Determining anchor residues in light and heavy chain by whether their index matches an IMGT-defined anchor-position
    for residue in model['L']:
        if residue.id[1] in IMGT_LC_ANCHORS:
            anchor_residues.add(residue)
    for residue in model['H']:
        if residue.id[1] in IMGT_HC_ANCHORS:
            anchor_residues.add(residue)
            
    return anchor_residues


def surface_exposed(model: PDB.Model.Model | PDB.Chain.Chain) -> set[PDB.Residue.Residue]:
    """
    Takes a PDB model or chain, returns a set of residues that have a sidechain relative surface exposure of at least 7.5%
    (relative to Ala-Arg-Ala).
    
    :param model: A PDB model or chain (invoked with ``.get_residues()``).
    :return surface_residues: A set of surface exposed residues.
    """
    # Computing Surface-Accessible Surface Area (SASA) using the Shrake-Rupley algorithm at atomic level
    sr = ShrakeRupley()
    sr.compute(model, level="A")
    
    # Defining a minimum relative exposed surface area threshold of 7.5%
    exposure_threshold = 0.075
    
    # Looping over residues
    surface_residues = set()
    for residue in model.get_residues():
        
        # Calculating SASA for sidechain atoms only #TODO: Relative SASA for entire residue or just side chain atoms?
        sasa = sum([atom.sasa for atom in residue if atom.name not in backbone_atoms])
        
        # Computing relative surface exposure relative to an Ala-Arg-Ala tripeptide
        relative_sasa = sasa / ARA_SASA
        
        # Determining whether a residue meets the surface exposure threshold
        if relative_sasa > exposure_threshold: 
            surface_residues.add(residue)
            
    return surface_residues
    

def cdr_vicinity(model: PDB.Model.Model) -> set[PDB.Residue.Residue]:
    """
    Takes a PDB model or chain, returns a set of surface-exposed residues that are
        1) CDRs,
        2) anchor residues, or
        3) in the vicinity of CDRs by a minimum distance of 4 Å between the closest heavy atoms.
    
    :param model: A PDB model.
    :return cdr_vicinity_residues: A set of surface exposed residues in, immediately flanking, or in the vicinity of CDRs.
    """
    # Initializing set with CDR and anchor residues
    cdr_vicinity_residues = cdr(model) | anchor(model)
    
    # Looping over surface exposed residues
    for res in surface_exposed(model):
        
        # Looping over possible neighbouring CDR residues
        for cdr_res in cdr(model):
            
            # Skipping to next surface-exposed residue if it is already assigned as in-vicinity
            if res in cdr_vicinity_residues:
                break
            
            # Looping over surface-exposed residue heavy atoms
            for atom in [atom for atom in res if atom.element != 'H']:
                
                # Skipping to next surface-exposed residue if it is already assigned as in-vicinity
                if res in cdr_vicinity_residues:
                    break
                
                # Looping over possible neighbouring CDR atoms
                for cdr_atom in [atom for atom in cdr_res if atom.element != 'H']:
                    
                    # Adding surface-exposed residue to set if in vicinity of a CDR residue
                    if atom - cdr_atom < 4.0:
                        cdr_vicinity_residues.add(res)
                        break
                    
    return cdr_vicinity_residues


def close_pairs(residues: set[PDB.Residue.Residue],) -> dict[tuple[PDB.Residue.Residue]: float]:
    """
    Takes a set of residues, returns a dictionary of distinct residue pair combinations and their minimum distance,
    for residue pairs with a minimum distance of at most 7.5 Å.
    
    :param residues: A set of PDB residues.
    :return close_residue_pairs: A dictionary with keys of tuples of pairs of residues and values of minimum distance between the residues.
    """
    # Initializing dictionary
    close_residue_pairs = {}
    
    # Looping over unique residue pair combinations
    for residue1, residue2 in itertools.combinations(residues, 2):
        
        # Keeping track of minimum distance
        min_dist = np.inf
        
        # Computing distance between each atom pair
        for atom1 in [atom for atom in residue1 if atom.element != 'H']:
            for atom2 in [atom for atom in residue2 if atom.element != 'H']:
                dist = atom1 - atom2
                
                # Noting new minimum distance
                if dist < min_dist:
                    min_dist = dist
                    
        # Determining whether residue pair meets minimum distance threshold
        if min_dist < 7.5:
            close_residue_pairs[(residue1, residue2)] = min_dist
        
    return close_residue_pairs


def salt_bridges(model: PDB.Model.Model) -> set[PDB.Residue.Residue]:
    """
    Takes a PDB model, returns a set of residues involved in salt bridges in the model.
    The definition for a salt bridge is for the negatively charged O of a Asp/Glu to be within 3.2 Å of a positively charged N of a Lys/Arg.
    
    :param model: A PDB model.
    :return salt_bridge_residues: A set of PDB residues that participate in salt bridges
    """
    # Initializing set
    salt_bridge_residues = set()
    
    # Identifying potential salt bridge residues
    positively_sb_residues = set([residue for residue in model.get_residues() if residue.resname in ['LYS', 'ARG']])
    negatively_sb_residues = set([residue for residue in model.get_residues() if residue.resname in ['ASP', 'GLU']])
    
    # Looping over all possible combinations
    for positive_residue, negative_residue in itertools.product(positively_sb_residues, negatively_sb_residues):
        
        # Identifying the positively charged N
        if positive_residue.resname == 'LYS':
            pos_atoms = [positive_residue['NZ']]
        elif positive_residue.resname == 'ARG':
            pos_atoms = [positive_residue['NE'], positive_residue['NH1'], positive_residue['NH2']] # Assuming any of the amines could be protonated
        
        # Identifying the negatively charged O
        if negative_residue.resname == 'ASP':
            neg_atoms = [negative_residue['OD1'], negative_residue['OD1']] # Either O could be negatively charged due to resonance
        elif negative_residue.resname == 'GLU':
            neg_atoms = [negative_residue['OE1'], negative_residue['OE1']] # Either O could be negatively charged due to resonance
            
        # Determining if the atomic distance is within threshold
        for pos_atom in pos_atoms:
            for neg_atom in neg_atoms:
                dist = pos_atom - neg_atom
                if dist <= 3.2:
                    salt_bridge_residues.add(positive_residue)
                    salt_bridge_residues.add(negative_residue)
        
        return salt_bridge_residues


# ················································································· #
    
def psh(residues: set[PDB.Residue.Residue]) -> float:
    """
    Takes a set of residues, computes the patches of surface hydrophobicity (PSH) score from Therapeutic Antibody Profiler (TAP).
    
    :param residues: A set of residues
    :return psh: The PSH score.
    """
    # Initializing sum
    psh = 0
    
    # Looping over close pairs of residues (i.e. with a minimum distance of 7.5 Å)
    close_residue_pairs = close_pairs(residues)
    for (residue1, residue2), dist in close_residue_pairs.items():
        
        # Retrievining normalized hydrophobicity scores
        H_1 = NORMALIZED_HYDROPHOBICITY.get(residue1.resname)
        H_2 = NORMALIZED_HYDROPHOBICITY.get(residue2.resname)
        
        # Adding pair value to PSH sum
        psh += (H_1 * H_2) / (dist**2)
        
    return psh


def ppc(residues: set[PDB.Residue.Residue], model: PDB.Model.Model) -> float:
    """
    Takes a set of residues, computes the patches of positive charge (PPC) score from Therapeutic Antibody Profiler (TAP).
    Does not take the charge contribution of salt bridge residues into account.
    
    :param residues: A set of residues
    :return ppc: The PPC score.
    """
    # Initializing sum
    ppc = 0
    
    # Skipping the contribution of salt bridge residues
    residues -= salt_bridges(model)
    
    # Looping over close pairs of residues (i.e. with a minimum distance of 7.5 Å)
    close_residue_pairs = close_pairs(residues)
    for (residue1, residue2), dist in close_residue_pairs.items():
        
        # Retrievining normalized hydrophobicity scores
        Q_1 = CHARGE.get(residue1.resname, 0)
        Q_2 = CHARGE.get(residue2.resname, 0)
        
        # Adding pair value to PPC sum if both residues are positively charged
        if Q_1 > 0 and Q_2 > 0:
            ppc += (Q_1 * Q_2) / (dist**2) # Always both positive, so analogous to absolute value
            
    return ppc


def pnc(residues: set[PDB.Residue.Residue], model: PDB.Model.Model) -> float:
    """
    Takes a set of residues, computes the patches of negative charge (PNC) score from Therapeutic Antibody Profiler (TAP).
    Does not take the charge contribution of salt bridge residues into account.
    
    :param residues: A set of residues
    :return pnc: The PNC score.
    """
    # Initializing sum
    pnc = 0
    
    # Skipping the contribution of salt bridge residues
    residues -= salt_bridges(model)
    
    # Looping over close pairs of residues (i.e. with a minimum distance of 7.5 Å)
    close_residue_pairs = close_pairs(residues)
    for (residue1, residue2), dist in close_residue_pairs.items():
        
        # Retrievining normalized hydrophobicity scores
        Q_1 = CHARGE.get(residue1.resname, 0)
        Q_2 = CHARGE.get(residue2.resname, 0)
        
        # Adding pair value to PPC sum if both residues are negative charged
        if Q_1 < 0 and Q_2 < 0:
            pnc += (Q_1 * Q_2) / (dist**2) # Always both negative, so analogous to absolute value
            
    return pnc


def sfvcsp(model: PDB.Model.Model) -> float:
    """
    Takes an PDB model with a heavy ``H`` and light ``L`` chain, returns the structural Fv charge symmetry parameter (SFvCSP) from Therapeutic Antibody Profiler (TAP).
    
    :param model: A PDB model with chains ``H`` and ``L``.
    :return sfvcsp: The SFvCSP score.
    """
    # Summing charge of surface-exposed residues for each chain
    sigmaLC_Q = sum([CHARGE.get(residue.resname, 0) for residue in surface_exposed(model['L'])])
    sigmaHC_Q = sum([CHARGE.get(residue.resname, 0) for residue in surface_exposed(model['H'])])
    
    # Computing SFvCSP measure
    sfvcsp = sigmaLC_Q * sigmaHC_Q
    
    return sfvcsp


# ················································································· #

def calc_tap(input_pdbs: list[str],
             output_tsv: str=None) -> pd.DataFrame | None:
    """
    Takes one or more paths to ABodyBuilder2-predicted antibody structures as PDB files,
    returns a table of [Therapeutic Antibody Profiler (TAP)](https://doi.org/10.1073/pnas.1810576116) measures from the Charlotte Dean lab.
    
    :param input_pdbs: A list of PDB paths. PDB files must have a ``H`` and ``L`` chain and have residue indexing that matches IMGT numbering.
    :param output_tsv: The path to write output to as TSV. If none is specified, a dataframe will be returned.
    :return: If ``output_tsv`` is not specified, a pandas dataframe of computed TAP measures.
    """
        
    # Looping over models
    tap = []
    for input_pdb in input_pdbs:
        
        # Loading ABodyBuilder2 model
        model = load_pdb(input_pdb)

        # Determining residues of interest
        cdr_residues = cdr(model)
        cdr_vicinity_residues = cdr_vicinity(model)
        
        # Calculating TAP measures
        tap.append(pd.Series({
            'cdr_len': len(cdr_residues),
            'psh_cdr': psh(cdr_vicinity_residues),
            'ppc_cdr': ppc(cdr_vicinity_residues, model),
            'pnc_cdr': pnc(cdr_vicinity_residues, model),
            'sfvcsp': sfvcsp(model),
        }, name=input_pdb.replace('.pdb', '')))
        
    # Combining results from all models
    tap = pd.DataFrame(tap)
    tap.index.name = 'model'
    
    # Returning or writing results
    if output_tsv is None:
        return tap
    else:
        tap.to_csv(output_tsv, sep='\t')


# ················································································· #

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='calc-tap',
        usage='''
            python calc-tap.py
                [-h]
                --input-pdbs         INPUT1.PDB [INPUT2.PDB ...]
            [   --output-tsv         OUTPUT.TSV  ]
                ''',
        description="""
        Calculate Therapeutic Antibody Profiler (TAP) measures, including
            total CDR length (cdr_len),
            patches of surface hydrophobicity in CDR vicinty (psh_cdr),
            patches of positive charge in CDR vicinty (ppc_cdr),
            patches of negative charge in CDR vicinty (pnc_cdr),
            and the structural Fv charge symmetry parameter (sfvcsp)
        from an ABodyBuilder2-generated PDB file.""",
        epilog='// Generated Script')

    parser.add_argument('--input-pdbs',
                        type=str,
                        nargs='+',
                        required=True,
                        help="Path(s) to ABodyBuilder2-generated PDB file(s).")
    parser.add_argument('--output-tsv',
                        type=str,
                        required=False,
                        default='tap.tsv',
                        help="The path to write the output TSV file to. (Defaults to tap.tsv)")
    
    # ············································································· #

    args = parser.parse_args()
    
    calc_tap(**vars(args))
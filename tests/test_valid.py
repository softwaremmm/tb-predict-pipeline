"""Suite of unit tests based on the `gnomonicus` test suite
"""
import json

import gnomonicus
import pytest

# Helpful function import for testing nested JSON equality as it gives exact differences
from recursive_diff import recursive_eq

"""
Due to complications testing equalities of nested jsons of lists/dicts, there is a lot of 
code specificially dedicated to ensuring the lists are in the same order (they differ due to
a mixture of dictionary behaviour and different positions within files). However, we only care that
the contents of the JSON is there for these tests rather than caring about order.
"""


def prep_json(j: dict) -> dict:
    """Prepare a JSON for comparison by removing fields which cannot be reproduced

    Args:
        j (dict): Initial JSON

    Returns:
        dict: JSON without fields such as time and file paths
    """
    del j["meta"]["time_taken_s"]
    del j["meta"]["UTC-datetime-completed"]
    del j["meta"]["catalogue_file"]
    del j["meta"]["reference_file"]
    del j["meta"]["vcf_file"]
    return j


def variants_key(x):
    """Used as the sorted(key=) function for reliably sorting the variants list

    Args:
        x (list): List of the ordered values as key, value pairs

    Returns:
        str: String of the `variant+gene_name`
    """
    variant = ""
    gene = ""
    for i in x:
        if i[0] == "variant":
            variant = i[1]
        elif i[0] == "gene_name":
            gene = i[1]
    return variant + gene


def ordered(obj):
    """Recursively sort a JSON for equality checking. Based on https://stackoverflow.com/questions/25851183/how-to-compare-two-json-objects-with-the-same-elements-in-a-different-order-equa

    Args:
        obj (object): Any JSON element. Probably one of dict, list, tuple, str, int, None

    Returns:
        object: Sorted JSON
    """
    if isinstance(obj, dict):
        if "variants" in obj.keys():
            # We have the 'data' field which needs a little extra nudge
            if "effects" in obj.keys():
                # Case when we have effects populated
                return [
                    ("antibiogram", ordered(obj["antibiogram"])),
                    ("effects", ordered(obj["effects"])),
                    ("mutations", ordered(obj["mutations"])),
                    (
                        "variants",
                        sorted([ordered(x) for x in obj["variants"]], key=variants_key),
                    ),
                ]
            elif "mutations" in obj.keys():
                # Case for if we have mutations and variants but no effects
                return [
                    ("mutations", ordered(obj["mutations"])),
                    (
                        "variants",
                        sorted([ordered(x) for x in obj["variants"]], key=variants_key),
                    ),
                ]
            else:
                # Case where no effects/mutations have been populated
                return [
                    (
                        "variants",
                        sorted([ordered(x) for x in obj["variants"]], key=variants_key),
                    )
                ]
        else:
            return sorted((k, ordered(obj[k])) for k in sorted(list(obj.keys())))

    if isinstance(obj, list) or isinstance(obj, tuple):
        return sorted(ordered(x) for x in obj)

    # Nones cause issues with ordering as there is no < operator.
    # Convert to string to avoid this
    if type(obj) == type(None):
        return str(obj)

    # Because nan types are helpful, `float('nan') == float('nan') -> False`
    # So check based on str value rather than direct equality
    if str(obj) == "nan":
        # Conversion to None for reproducability
        return str(None)

    if isinstance(obj, int):
        # Ints are still ordered (just not numerically) if sorted by str value, so convert to str
        # to allow sorting lists of None/int
        return str(obj)
    if isinstance(obj, float):
        # Similarly convert float, but check if they are x.0 to avoid str comparison issues
        if int(obj) == obj:
            return str(int(obj))
        else:
            return str(obj)
    else:
        return obj


def test_1():
    """Input:
        NC_045512.2-S_E484K-minos.vcf
    Expect output:
        variants:    23012g>a
        mutations:   S@E484K
        predictions: {'AAA': 'R', 'BBB': 'S'}
    """
    vcfStem = "NC_045512.2-S_E484K-minos"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "REF": "g",
                        "ALTS": ["a"],
                        "POS": 23012,
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [
                {
                    "mutation": "E484K",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "aaa",
                }
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "E484K",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/1/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_3():
    """Input:
        NC_045512.2-S_F2F-minos.vcf
    Expect output:
        variants:    21568t>c
        mutations:   S@F2F
        predictions: {'AAA': 'S', 'BBB': 'S'}
    """
    vcfStem = "NC_045512.2-S_F2F-minos"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "21568t>c",
                    "nucleotide_index": 21568,
                    "gene_name": "S",
                    "gene_position": 2,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21568,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [
                {
                    "mutation": "F2F",
                    "gene": "S",
                    "gene_position": 2,
                    "ref": "ttt",
                    "alt": "ttc",
                },
                {"mutation": "t6c", "gene": "S", "gene_position": 6},
            ],
            "effects": {
                "AAA": [
                    {"gene": "S", "mutation": "F2F", "prediction": "S", "evidence": {}},
                    {"phenotype": "S"},
                ],
            },
            "antibiogram": {"AAA": "S", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/3/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_4():
    """Input:
        NC_045512.2-S_F2L-minos.vcf
    Expect output:
        variants:    21566t>c
        mutations:   S@F2L
        predictions: {'AAA': 'U', 'BBB': 'S'}
    """
    # Setup
    vcfStem = "NC_045512.2-S_F2L-minos"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "21566t>c",
                    "nucleotide_index": 21566,
                    "gene_name": "S",
                    "gene_position": 2,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21566,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [
                {
                    "mutation": "F2L",
                    "gene": "S",
                    "gene_position": 2,
                    "ref": "ttt",
                    "alt": "ctt",
                }
            ],
            "effects": {
                "AAA": [
                    {"gene": "S", "mutation": "F2L", "prediction": "U", "evidence": {}},
                    {"phenotype": "U"},
                ],
            },
            "antibiogram": {"AAA": "U", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/4/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_5():
    """Input:
        NC_045512.2-S_200_indel-minos.vcf
    Expect output:
        variants:    21762_ins_c
        mutations:   S@200_ins_c
        predictions: {'AAA': 'R', 'BBB': 'S'}
    """
    vcfStem = "NC_045512.2-S_200_indel-minos"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "21762_ins_c",
                    "nucleotide_index": 21762,
                    "gene_name": "S",
                    "gene_position": 200,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21762,
                        "REF": "c",
                        "ALTS": ["cc"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "200_ins_c", "gene": "S", "gene_position": 200},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "200_ins_c",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/5/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_6():
    """Input:
        NC_045512.2-double-minos.vcf
    Expect output:
        variants:    27758g>c
        mutations:   ORF7a!122S, ORF7b@M1I
        predictions: {'AAA': 'R', 'BBB': 'R'}
    """
    vcfStem = "NC_045512.2-double-minos"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "27758g>c",
                    "nucleotide_index": 27758,
                    "gene_name": "ORF7a",
                    "gene_position": 122,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "27758g>c",
                    "nucleotide_index": 27758,
                    "gene_name": "ORF7b",
                    "gene_position": 1,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "!122S",
                    "gene": "ORF7a",
                    "gene_position": 122,
                    "ref": "tga",
                    "alt": "tca",
                },
                {
                    "mutation": "M1I",
                    "gene": "ORF7b",
                    "gene_position": 1,
                    "ref": "atg",
                    "alt": "atc",
                },
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "ORF7a",
                        "mutation": "!122S",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "BBB": [
                    {
                        "gene": "ORF7b",
                        "mutation": "M1I",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "R"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/6/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_7():
    """Input:
        NC_045512.2-S_E484K&1450_ins_a-minos.vcf
    Expect output:
        variants:    23012g>a, 23012_ins_a
        mutations:   S@E484K, S@1450_ins_a, S@1450_ins_a&S@E484K
        predictions: {'AAA': 'R', 'BBB': 'R'}
    """
    vcfStem = "NC_045512.2-S_E484K&1450_ins_a-minos"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "23012_ins_a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 1450,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["aa"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["aa"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "1450_ins_a", "gene": "S", "gene_position": 1450},
                {
                    "mutation": "E484K",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "aaa",
                },
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "E484K",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "1450_ins_a",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "BBB": [
                    {
                        "gene": None,
                        "mutation": "S@1450_ins_a&S@E484K",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "R"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/7/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_8():
    """Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:0.068, 25283_del_g:0.068, 25252_ins_cc:0.068
        mutations:   !1274Q:0.068, 3721_del_g:0.068, 3690_ins_cc:0.068
        predictions: {'AAA': 'R'}
    """
    vcfStem = "NC_045512.2-minors"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "25382t>c:0.067",
                    "nucleotide_index": 25382,
                    "gene_name": "S",
                    "gene_position": 1274,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25382,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "21558g>a:0.067",
                    "nucleotide_index": 21558,
                    "gene_name": "S",
                    "gene_position": -5,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21558,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25252_ins_cc:0.067",
                    "nucleotide_index": 25252,
                    "gene_name": "S",
                    "gene_position": 3690,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25252,
                        "REF": "g",
                        "ALTS": ["gcc"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25283_del_t:0.067",
                    "nucleotide_index": 25283,
                    "gene_name": "S",
                    "gene_position": 3721,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25282,
                        "REF": "tt",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "!1274Q:0.067",
                    "gene": "S",
                    "gene_position": 1274,
                    "ref": "taa",
                    "alt": "caa",
                },
                {"mutation": "g-5a:0.067", "gene": "S", "gene_position": -5},
                {"mutation": "3721_del_t:0.067", "gene": "S", "gene_position": 3721},
                {"mutation": "3690_ins_cc:0.067", "gene": "S", "gene_position": 3690},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "!1274Q:0.067",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "g-5a:0.067",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "3721_del_t:0.067",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "3690_ins_cc:0.067",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/8/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_9():
    """Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:3, 25283_del_g:3, 25252_ins_cc:3
        mutations:   !1274Q:3, 3721_del_g:3, 3690_ins_cc:3
        predictions: {'AAA': 'R'}
    """
    vcfStem = "NC_045512.2-minors"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "25382t>c:3",
                    "nucleotide_index": 25382,
                    "gene_name": "S",
                    "gene_position": 1274,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25382,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "21558g>a:3",
                    "nucleotide_index": 21558,
                    "gene_name": "S",
                    "gene_position": -5,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21558,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25252_ins_cc:3",
                    "nucleotide_index": 25252,
                    "gene_name": "S",
                    "gene_position": 3690,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25252,
                        "REF": "g",
                        "ALTS": ["gcc"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25283_del_t:3",
                    "nucleotide_index": 25283,
                    "gene_name": "S",
                    "gene_position": 3721,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 3],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25282,
                        "REF": "tt",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "!1274Q:3",
                    "gene": "S",
                    "gene_position": 1274,
                    "ref": "taa",
                    "alt": "caa",
                },
                {"mutation": "g-5a:3", "gene": "S", "gene_position": -5},
                {"mutation": "3721_del_t:3", "gene": "S", "gene_position": 3721},
                {"mutation": "3690_ins_cc:3", "gene": "S", "gene_position": 3690},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "!1274Q:3",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "g-5a:3",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "3721_del_t:3",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "3690_ins_cc:3",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/9/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_10():
    """Testing a catalogue and sample which have large deletions
    Input:
        TEST-DNA-large-del.vcf
    Expect output:
        variants:    3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc
        mutations:   A@-1_del_aaaaaaaaccccccccccgggggggggg, A@del_0.93, B@del_1.0, C@4_del_ggg
        predictions: {'AAA': 'R'}
    """
    vcfStem = "TEST-DNA-large-del"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "TEST_DNA",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test_dna",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                    "nucleotide_index": 3,
                    "gene_name": "A",
                    "gene_position": -1,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 4,
                        "COV": [1, 3],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                    "nucleotide_index": 3,
                    "gene_name": "B",
                    "gene_position": 1,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 4,
                        "COV": [1, 3],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                    "nucleotide_index": 3,
                    "gene_name": "C",
                    "gene_position": 4,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 4,
                        "COV": [1, 3],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "-1_del_aaaaaaaaccccccccccgggggggggg",
                    "gene": "A",
                    "gene_position": -1,
                },
                {"mutation": "del_0.93", "gene": "A", "gene_position": None},
                {"mutation": "del_1.0", "gene": "B", "gene_position": None},
                {
                    "gene": "B",
                    "gene_position": 1,
                    "mutation": "1_del_gggttttttttttaaaaaaaaaacccccccccc",
                },
                {"mutation": "4_del_ggg", "gene": "C", "gene_position": 4},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "A",
                        "mutation": "-1_del_aaaaaaaaccccccccccgggggggggg",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "A",
                        "mutation": "del_0.93",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "B",
                        "mutation": "del_1.0",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "B",
                        "mutation": "1_del_gggttttttttttaaaaaaaaaacccccccccc",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "C",
                        "mutation": "4_del_ggg",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {
                "AAA": "R",
            },
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/10/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


# ===================================================
#                   TB TESTS
# ===================================================


def test_11():
    vcfStem = "NC_000962_3_test_0001"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_000962",
            "catalogue_type": "RFUS",
            "catalogue_name": "test_001",
            "catalogue_version": "v1.00",
        },
        "data": {
            "variants": [
                {
                    "variant": "7571g>a",
                    "nucleotide_index": 7571,
                    "gene_name": "gyrA",
                    "gene_position": 90,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 54.0,
                        "ALLELE_DP": [0.0, 54.0],
                        "FRS": 1.0,
                        "COV_TOTAL": 54,
                        "COV": [0, 54],
                        "GT_CONF": 31.81,
                        "GT_CONF_PERCENTILE": 11.47,
                        "POS": 7571,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "7585g>c",
                    "nucleotide_index": 7585,
                    "gene_name": "gyrA",
                    "gene_position": 95,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 54.0,
                        "ALLELE_DP": [0.0, 54.0],
                        "FRS": 1.0,
                        "COV_TOTAL": 54,
                        "COV": [0, 54],
                        "GT_CONF": 31.81,
                        "GT_CONF_PERCENTILE": 11.47,
                        "POS": 7585,
                        "REF": "g",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "760854a>t",
                    "nucleotide_index": 760854,
                    "gene_name": "rpoB",
                    "gene_position": 350,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 43.0,
                        "ALLELE_DP": [0.0, 43.0],
                        "FRS": 1.0,
                        "COV_TOTAL": 43,
                        "COV": [0, 43],
                        "GT_CONF": 396.18,
                        "GT_CONF_PERCENTILE": 38.85,
                        "POS": 760854,
                        "REF": "a",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "761155c>t",
                    "nucleotide_index": 761155,
                    "gene_name": "rpoB",
                    "gene_position": 450,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 761155,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1473183a>c",
                    "nucleotide_index": 1473183,
                    "gene_name": "rrs",
                    "gene_position": 1338,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 209.84,
                        "GT_CONF_PERCENTILE": 77.55,
                        "POS": 1473183,
                        "REF": "a",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1473246a>g",
                    "nucleotide_index": 1473246,
                    "gene_name": "rrs",
                    "gene_position": 1401,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 43.89,
                        "GT_CONF_PERCENTILE": 23.15,
                        "POS": 1473246,
                        "REF": "a",
                        "ALTS": ["g"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1674048g>t",
                    "nucleotide_index": 1674048,
                    "gene_name": "fabG1",
                    "gene_position": 203,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 54.0,
                        "ALLELE_DP": [4.0, 50.0],
                        "FRS": 0.9259,
                        "COV_TOTAL": 54,
                        "COV": [4, 50],
                        "GT_CONF": 98.34,
                        "GT_CONF_PERCENTILE": 99.05,
                        "POS": 1674048,
                        "REF": "g",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289010c>a",
                    "nucleotide_index": 2289010,
                    "gene_name": "pncA",
                    "gene_position": 78,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "POS": 2289010,
                        "REF": "c",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289193c>t",
                    "nucleotide_index": 2289193,
                    "gene_name": "pncA",
                    "gene_position": 17,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 0.57,
                        "GT_CONF_PERCENTILE": 85.13,
                        "POS": 2289193,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "a1338c", "gene": "rrs", "gene_position": 1338},
                {"mutation": "a1401g", "gene": "rrs", "gene_position": 1401},
                {
                    "mutation": "T350S",
                    "gene": "rpoB",
                    "gene_position": 350,
                    "ref": "acc",
                    "alt": "tcc",
                },
                {
                    "mutation": "S450L",
                    "gene": "rpoB",
                    "gene_position": 450,
                    "ref": "tcg",
                    "alt": "ttg",
                },
                {
                    "mutation": "A90A",
                    "gene": "gyrA",
                    "gene_position": 90,
                    "ref": "gcg",
                    "alt": "gca",
                },
                {"mutation": "g270a", "gene": "gyrA", "gene_position": 270},
                {
                    "mutation": "S95T",
                    "gene": "gyrA",
                    "gene_position": 95,
                    "ref": "agc",
                    "alt": "acc",
                },
                {
                    "mutation": "G17S",
                    "gene": "pncA",
                    "gene_position": 17,
                    "ref": "ggc",
                    "alt": "agc",
                },
                {
                    "mutation": "G78C",
                    "gene": "pncA",
                    "gene_position": 78,
                    "ref": "ggc",
                    "alt": "tgc",
                },
                {
                    "mutation": "L203L",
                    "gene": "fabG1",
                    "gene_position": 203,
                    "ref": "ctg",
                    "alt": "ctt",
                },
                {"mutation": "g609t", "gene": "fabG1", "gene_position": 609},
            ],
            "effects": {
                "KAN": [
                    {
                        "gene": "rrs",
                        "mutation": "a1338c",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "rrs",
                        "mutation": "a1401g",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "RIF": [
                    {
                        "gene": "rpoB",
                        "mutation": "T350S",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "rpoB",
                        "mutation": "S450L",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "MXF": [
                    {
                        "gene": "gyrA",
                        "mutation": "A90A",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "gyrA",
                        "mutation": "S95T",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "S"},
                ],
                "PZA": [
                    {
                        "gene": "pncA",
                        "mutation": "G17S",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "G78C",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "INH": [
                    {
                        "gene": "fabG1",
                        "mutation": "L203L",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"KAN": "R", "RIF": "R", "MXF": "S", "PZA": "R", "INH": "R"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/11/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_12():
    vcfStem = "NC_000962_3_test_0002"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_000962",
            "catalogue_type": "RFUS",
            "catalogue_name": "test_001",
            "catalogue_version": "v1.00",
        },
        "data": {
            "antibiogram": {"INH": "R", "KAN": "F", "MXF": "U", "PZA": "R", "RIF": "U"},
            "variants": [
                {
                    "variant": "1674048g>x",
                    "nucleotide_index": 1674048,
                    "gene_name": "fabG1",
                    "gene_position": 203,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "COV": [0, 0],
                        "COV_TOTAL": 0,
                        "GT_CONF_PERCENTILE": 0.0,
                        "GT": [None, None],
                        "FRS": ".",
                        "DP": None,
                        "ALLELE_DP": [0, 0],
                        "GT_CONF": 0.0,
                        "POS": 1674048,
                        "REF": "g",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": None,
                },
                {
                    "variant": "7570c>t:0.05",
                    "nucleotide_index": 7570,
                    "gene_name": "gyrA",
                    "gene_position": 90,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "FRS": 1.0,
                        "GT": [0, 0],
                        "ALLELE_DP": [95, 5],
                        "COV_TOTAL": 100,
                        "DP": 100,
                        "COV": [95, 5],
                        "GT_CONF": 31.81,
                        "GT_CONF_PERCENTILE": 11.47,
                        "POS": 7570,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "7581g>a",
                    "nucleotide_index": 7581,
                    "gene_name": "gyrA",
                    "gene_position": 94,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 54,
                        "ALLELE_DP": [0, 54],
                        "GT_CONF": 31.81,
                        "GT_CONF_PERCENTILE": 11.47,
                        "COV_TOTAL": 54,
                        "FRS": 1.0,
                        "COV": [0, 54],
                        "POS": 7581,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "7585g>c",
                    "nucleotide_index": 7585,
                    "gene_name": "gyrA",
                    "gene_position": 95,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "DP": 50,
                        "FRS": 0.98,
                        "GT_CONF_PERCENTILE": 11.47,
                        "COV": [0, 1, 49],
                        "GT_CONF": 31.81,
                        "COV_TOTAL": 100,
                        "ALLELE_DP": [0, 1, 49],
                        "GT": [2, 2],
                        "POS": 7585,
                        "REF": "g",
                        "ALTS": ["a", "c"],
                    },
                    "vcf_idx": 2,
                },
                {
                    "variant": "2155168c>g",
                    "nucleotide_index": 2155168,
                    "gene_name": "katG",
                    "gene_position": 315,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "DP": 89,
                        "GT_CONF": 204.17,
                        "ALLELE_DP": [5, 84],
                        "GT": [1, 1],
                        "FRS": 0.9438,
                        "COV": [5, 84],
                        "GT_CONF_PERCENTILE": 49.86,
                        "COV_TOTAL": 89,
                        "POS": 2155168,
                        "REF": "c",
                        "ALTS": ["g"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2154543c>g",
                    "nucleotide_index": 2154543,
                    "gene_name": "katG",
                    "gene_position": 523,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "DP": 89,
                        "COV_TOTAL": 89,
                        "GT": [1, 1],
                        "FRS": 0.9438,
                        "ALLELE_DP": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "COV": [5, 84],
                        "POS": 2154543,
                        "REF": "c",
                        "ALTS": ["g"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289252t>c",
                    "nucleotide_index": 2289252,
                    "gene_name": "pncA",
                    "gene_position": -11,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "ALLELE_DP": [5, 84],
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT": [1, 1],
                        "DP": 89,
                        "FRS": 0.9438,
                        "POS": 2289252,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289237_del_cg",
                    "nucleotide_index": 2289237,
                    "gene_name": "pncA",
                    "gene_position": 4,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT_CONF": 204.17,
                        "GT": [1, 1],
                        "DP": 89,
                        "ALLELE_DP": [5, 84],
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF_PERCENTILE": 49.86,
                        "FRS": 0.9438,
                        "POS": 2289236,
                        "REF": "ccg",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2288945g>t",
                    "nucleotide_index": 2288945,
                    "gene_name": "pncA",
                    "gene_position": 99,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "FRS": 0.9438,
                        "DP": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "COV_TOTAL": 89,
                        "ALLELE_DP": [5, 84],
                        "GT": [1, 1],
                        "POS": 2288945,
                        "REF": "g",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2288847_ins_ttt",
                    "nucleotide_index": 2288847,
                    "gene_name": "pncA",
                    "gene_position": 394,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "DP": 89,
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "GT_CONF": 0.57,
                        "GT": [1, 1],
                        "COV": [5, 84],
                        "GT_CONF_PERCENTILE": 85.13,
                        "ALLELE_DP": [5, 84],
                        "POS": 2288847,
                        "REF": "c",
                        "ALTS": ["cttt"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "760854a>c",
                    "nucleotide_index": 760854,
                    "gene_name": "rpoB",
                    "gene_position": 350,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "ALLELE_DP": [0, 43],
                        "GT": [1, 1],
                        "COV_TOTAL": 43,
                        "DP": 43,
                        "GT_CONF": 396.18,
                        "COV": [0, 43],
                        "FRS": 1.0,
                        "GT_CONF_PERCENTILE": 38.85,
                        "POS": 760854,
                        "REF": "a",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "761156g>a",
                    "nucleotide_index": 761156,
                    "gene_name": "rpoB",
                    "gene_position": 450,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "ALLELE_DP": [2, 98],
                        "COV_TOTAL": 100,
                        "GT": [1, 1],
                        "DP": 100,
                        "FRS": 0.98,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 761156,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "761160c>x",
                    "nucleotide_index": 761160,
                    "gene_name": "rpoB",
                    "gene_position": 452,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "COV": [0, 0],
                        "COV_TOTAL": 0,
                        "GT_CONF_PERCENTILE": 0.0,
                        "GT_CONF": 0.0,
                        "DP": None,
                        "ALLELE_DP": [0, 0],
                        "GT": [None, None],
                        "FRS": ".",
                        "POS": 761160,
                        "REF": "c",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": None,
                },
                {
                    "variant": "1473183a>x",
                    "nucleotide_index": 1473183,
                    "gene_name": "rrs",
                    "gene_position": 1338,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "DP": None,
                        "COV": [0, 0],
                        "FRS": ".",
                        "ALLELE_DP": [0, 0],
                        "GT_CONF_PERCENTILE": 0.0,
                        "COV_TOTAL": 0,
                        "GT_CONF": 0.0,
                        "GT": [None, None],
                        "POS": 1473183,
                        "REF": "a",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": None,
                },
                {
                    "variant": "1473245c>t",
                    "nucleotide_index": 1473245,
                    "gene_name": "rrs",
                    "gene_position": 1400,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100,
                        "ALLELE_DP": [2, 98],
                        "GT_CONF": 43.89,
                        "COV_TOTAL": 100,
                        "FRS": 0.98,
                        "GT_CONF_PERCENTILE": 23.15,
                        "COV": [2, 98],
                        "POS": 1473245,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "gene": "fabG1",
                    "gene_position": 203,
                    "mutation": "L203X",
                    "ref": "ctg",
                    "alt": "ctx",
                },
                {
                    "gene": "gyrA",
                    "gene_position": 90,
                    "mutation": "A90V:0.05",
                    "ref": "gcg",
                    "alt": "gtg",
                },
                {
                    "gene": "gyrA",
                    "gene_position": 94,
                    "mutation": "D94N",
                    "ref": "gac",
                    "alt": "aac",
                },
                {
                    "gene": "gyrA",
                    "gene_position": 95,
                    "mutation": "S95T",
                    "ref": "agc",
                    "alt": "acc",
                },
                {
                    "gene": "katG",
                    "gene_position": 315,
                    "mutation": "S315T",
                    "ref": "agc",
                    "alt": "acc",
                },
                {
                    "gene": "katG",
                    "gene_position": 523,
                    "mutation": "E523D",
                    "ref": "gag",
                    "alt": "gac",
                },
                {"gene": "pncA", "gene_position": -11, "mutation": "a-11g"},
                {"gene": "pncA", "gene_position": 4, "mutation": "4_del_cg"},
                {
                    "gene": "pncA",
                    "gene_position": 99,
                    "mutation": "Y99!",
                    "ref": "tac",
                    "alt": "taa",
                },
                {"gene": "pncA", "gene_position": 394, "mutation": "394_ins_aaa"},
                {
                    "gene": "rpoB",
                    "gene_position": 350,
                    "mutation": "T350P",
                    "ref": "acc",
                    "alt": "ccc",
                },
                {
                    "gene": "rpoB",
                    "gene_position": 450,
                    "mutation": "S450S",
                    "ref": "tcg",
                    "alt": "tca",
                },
                {
                    "gene": "rpoB",
                    "gene_position": 452,
                    "mutation": "L452X",
                    "ref": "ctg",
                    "alt": "xtg",
                },
                {"gene": "rpoB", "gene_position": 1350, "mutation": "g1350a"},
                {"gene": "rrs", "gene_position": 1338, "mutation": "a1338x"},
                {"gene": "rrs", "gene_position": 1400, "mutation": "c1400t"},
            ],
            "effects": {
                "INH": [
                    {
                        "gene": "fabG1",
                        "mutation": "L203X",
                        "prediction": "F",
                        "evidence": {},
                    },
                    {
                        "gene": "katG",
                        "mutation": "E523D",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "katG",
                        "mutation": "S315T",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "KAN": [
                    {
                        "gene": "rrs",
                        "mutation": "a1338x",
                        "prediction": "F",
                        "evidence": {},
                    },
                    {
                        "gene": "rrs",
                        "mutation": "c1400t",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {"phenotype": "F"},
                ],
                "MXF": [
                    {
                        "gene": "gyrA",
                        "mutation": "A90V:0.05",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "gyrA",
                        "mutation": "D94N",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "gyrA",
                        "mutation": "S95T",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "U"},
                ],
                "PZA": [
                    {
                        "gene": "pncA",
                        "mutation": "394_ins_aaa",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "4_del_cg",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "Y99!",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "a-11g",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "RIF": [
                    {
                        "gene": "rpoB",
                        "mutation": "S450S",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "rpoB",
                        "mutation": "T350P",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {"phenotype": "U"},
                ],
            },
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/12/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_13():
    vcfStem = "NC_000962_3_test_0003"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_000962",
            "catalogue_type": "RFUS",
            "catalogue_name": "test_001",
            "catalogue_version": "v1.00",
        },
        "data": {
            "variants": [
                {
                    "variant": "761154t>g",
                    "nucleotide_index": 761154,
                    "gene_name": "rpoB",
                    "gene_position": 450,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 761154,
                        "REF": "t",
                        "ALTS": ["g"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "761155c>t",
                    "nucleotide_index": 761155,
                    "gene_name": "rpoB",
                    "gene_position": 450,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 761155,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1673425c>t",
                    "nucleotide_index": 1673425,
                    "gene_name": "fabG1",
                    "gene_position": -15,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 43.89,
                        "GT_CONF_PERCENTILE": 23.15,
                        "POS": 1673425,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1673432t>a",
                    "nucleotide_index": 1673432,
                    "gene_name": "fabG1",
                    "gene_position": -8,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 43.89,
                        "GT_CONF_PERCENTILE": 23.15,
                        "POS": 1673432,
                        "REF": "t",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1674046c>g",
                    "nucleotide_index": 1674046,
                    "gene_name": "fabG1",
                    "gene_position": 203,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 43.89,
                        "GT_CONF_PERCENTILE": 23.15,
                        "POS": 1674046,
                        "REF": "c",
                        "ALTS": ["g"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289252t>a",
                    "nucleotide_index": 2289252,
                    "gene_name": "pncA",
                    "gene_position": -11,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "POS": 2289252,
                        "REF": "t",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "7402_del_gc",
                    "nucleotide_index": 7402,
                    "gene_name": "gyrA",
                    "gene_position": 101,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 7401,
                        "REF": "agc",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "761006_ins_tcgt",
                    "nucleotide_index": 761006,
                    "gene_name": "rpoB",
                    "gene_position": 1200,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 761006,
                        "REF": "c",
                        "ALTS": ["ctcgt"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1473183_ins_tcg",
                    "nucleotide_index": 1473183,
                    "gene_name": "rrs",
                    "gene_position": 1338,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 1473183,
                        "REF": "a",
                        "ALTS": ["atcg"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2156013_del_gc",
                    "nucleotide_index": 2156013,
                    "gene_name": "katG",
                    "gene_position": 98,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "POS": 2156012,
                        "REF": "cgc",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2288855_ins_tcg",
                    "nucleotide_index": 2288855,
                    "gene_name": "pncA",
                    "gene_position": 386,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "POS": 2288855,
                        "REF": "t",
                        "ALTS": ["ttcg"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289212_del_c",
                    "nucleotide_index": 2289212,
                    "gene_name": "pncA",
                    "gene_position": 30,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "POS": 2289211,
                        "REF": "tc",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289237_ins_aa",
                    "nucleotide_index": 2289237,
                    "gene_name": "pncA",
                    "gene_position": 4,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "POS": 2289237,
                        "REF": "c",
                        "ALTS": ["caa"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "7570c>t:0.1",
                    "nucleotide_index": 7570,
                    "gene_name": "gyrA",
                    "gene_position": 90,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 100.0,
                        "ALLELE_DP": [90.0, 10.0],
                        "FRS": 1.0,
                        "COV_TOTAL": 100,
                        "COV": [90, 10],
                        "GT_CONF": 31.81,
                        "GT_CONF_PERCENTILE": 11.47,
                        "POS": 7570,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "1338_ins_tcg", "gene": "rrs", "gene_position": 1338.0},
                {"mutation": "98_del_gc", "gene": "katG", "gene_position": 98.0},
                {"mutation": "a-11t", "gene": "pncA", "gene_position": -11.0},
                {"mutation": "4_ins_tt", "gene": "pncA", "gene_position": 4.0},
                {"mutation": "30_del_g", "gene": "pncA", "gene_position": 30.0},
                {"mutation": "386_ins_cga", "gene": "pncA", "gene_position": 386.0},
                {"mutation": "101_del_gc", "gene": "gyrA", "gene_position": 101.0},
                {
                    "mutation": "L203V",
                    "gene": "fabG1",
                    "gene_position": 203.0,
                    "ref": "ctg",
                    "alt": "gtg",
                },
                {"mutation": "c-15t", "gene": "fabG1", "gene_position": -15.0},
                {"mutation": "t-8a", "gene": "fabG1", "gene_position": -8.0},
                {
                    "mutation": "S450V",
                    "gene": "rpoB",
                    "gene_position": 450.0,
                    "ref": "tcg",
                    "alt": "gtg",
                },
                {"mutation": "1200_ins_tcgt", "gene": "rpoB", "gene_position": 1200.0},
                {
                    "mutation": "A90V:0.1",
                    "gene": "gyrA",
                    "gene_position": 90.0,
                    "ref": "gcg",
                    "alt": "gtg",
                },
            ],
            "effects": {
                "KAN": [
                    {
                        "gene": "rrs",
                        "mutation": "1338_ins_tcg",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {"phenotype": "U"},
                ],
                "INH": [
                    {
                        "gene": "katG",
                        "mutation": "98_del_gc",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "fabG1",
                        "mutation": "L203V",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "fabG1",
                        "mutation": "c-15t",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {
                        "gene": "fabG1",
                        "mutation": "t-8a",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "PZA": [
                    {
                        "gene": "pncA",
                        "mutation": "a-11t",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "4_ins_tt",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "30_del_g",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "386_ins_cga",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "MXF": [
                    {
                        "gene": "gyrA",
                        "mutation": "101_del_gc",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "gyrA",
                        "mutation": "A90V:0.1",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "RIF": [
                    {
                        "gene": "rpoB",
                        "mutation": "S450V",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "rpoB",
                        "mutation": "1200_ins_tcgt",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {"phenotype": "U"},
                ],
            },
            "antibiogram": {"KAN": "U", "INH": "R", "PZA": "R", "MXF": "R", "RIF": "U"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/13/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_14():
    """Input:
        NC_045512.2-S_E484K-minos.vcf with all N fasta
    Expect output:
        variants:    23012g>x
        mutations:   S@E484X
        predictions: {'AAA': 'F', 'BBB': 'S'}
    """
    vcfStem = "NC_045512.2-S_E484K-minos"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "REF": "g",
                        "ALTS": ["a"],
                        "POS": 23012,
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "23013a>x",
                    "nucleotide_index": 23013,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [None, None],
                        "COV": 0,
                        "POS": 23013,
                        "REF": "c",
                        "ALTS": [],
                    },
                    "vcf_idx": None,
                },
                {
                    "variant": "23014a>x",
                    "nucleotide_index": 23014,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [None, None],
                        "COV": 0,
                        "POS": 23014,
                        "REF": "c",
                        "ALTS": [],
                    },
                    "vcf_idx": None,
                },
            ],
            "mutations": [
                {
                    "mutation": "E484X",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "axx",
                },
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "E484X",
                        "prediction": "F",
                        "evidence": {},
                    },
                    {"phenotype": "F"},
                ],
            },
            "antibiogram": {"AAA": "F", "BBB": "S"},
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/14/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_15():
    # This should exactly match the outputs of test 11, but with additional `F` calls
    vcfStem = "NC_000962_3_test_0001"

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_000962",
            "catalogue_type": "RFUS",
            "catalogue_name": "test_001",
            "catalogue_version": "v1.00",
        },
        "data": {
            "variants": [
                {
                    "variant": "7571g>a",
                    "nucleotide_index": 7571,
                    "gene_name": "gyrA",
                    "gene_position": 90,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 54.0,
                        "ALLELE_DP": [0.0, 54.0],
                        "FRS": 1.0,
                        "COV_TOTAL": 54,
                        "COV": [0, 54],
                        "GT_CONF": 31.81,
                        "GT_CONF_PERCENTILE": 11.47,
                        "POS": 7571,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "7585g>c",
                    "nucleotide_index": 7585,
                    "gene_name": "gyrA",
                    "gene_position": 95,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 54.0,
                        "ALLELE_DP": [0.0, 54.0],
                        "FRS": 1.0,
                        "COV_TOTAL": 54,
                        "COV": [0, 54],
                        "GT_CONF": 31.81,
                        "GT_CONF_PERCENTILE": 11.47,
                        "POS": 7585,
                        "REF": "g",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "760854a>t",
                    "nucleotide_index": 760854,
                    "gene_name": "rpoB",
                    "gene_position": 350,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 43.0,
                        "ALLELE_DP": [0.0, 43.0],
                        "FRS": 1.0,
                        "COV_TOTAL": 43,
                        "COV": [0, 43],
                        "GT_CONF": 396.18,
                        "GT_CONF_PERCENTILE": 38.85,
                        "POS": 760854,
                        "REF": "a",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "761155c>t",
                    "nucleotide_index": 761155,
                    "gene_name": "rpoB",
                    "gene_position": 450,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 252.06,
                        "GT_CONF_PERCENTILE": 43.15,
                        "POS": 761155,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1473183a>c",
                    "nucleotide_index": 1473183,
                    "gene_name": "rrs",
                    "gene_position": 1338,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 209.84,
                        "GT_CONF_PERCENTILE": 77.55,
                        "POS": 1473183,
                        "REF": "a",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1473246a>g",
                    "nucleotide_index": 1473246,
                    "gene_name": "rrs",
                    "gene_position": 1401,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 100.0,
                        "ALLELE_DP": [2.0, 98.0],
                        "FRS": 0.98,
                        "COV_TOTAL": 100,
                        "COV": [2, 98],
                        "GT_CONF": 43.89,
                        "GT_CONF_PERCENTILE": 23.15,
                        "POS": 1473246,
                        "REF": "a",
                        "ALTS": ["g"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "1674046c>x",
                    "nucleotide_index": 1674046,
                    "gene_name": "fabG1",
                    "gene_position": 203,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [None, None],
                        "COV": 0,
                        "POS": 1674046,
                        "REF": "c",
                        "ALTS": [],
                    },
                    "vcf_idx": None,
                },
                {
                    "variant": "1674047t>x",
                    "nucleotide_index": 1674047,
                    "gene_name": "fabG1",
                    "gene_position": 203,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [None, None],
                        "COV": 0,
                        "POS": 1674047,
                        "REF": "c",
                        "ALTS": [],
                    },
                    "vcf_idx": None,
                },
                {
                    "variant": "1674048g>t",
                    "nucleotide_index": 1674048,
                    "gene_name": "fabG1",
                    "gene_position": 203,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 54.0,
                        "ALLELE_DP": [4.0, 50.0],
                        "FRS": 0.9259,
                        "COV_TOTAL": 54,
                        "COV": [4, 50],
                        "GT_CONF": 98.34,
                        "GT_CONF_PERCENTILE": 99.05,
                        "POS": 1674048,
                        "REF": "g",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289010c>a",
                    "nucleotide_index": 2289010,
                    "gene_name": "pncA",
                    "gene_position": 78,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 204.17,
                        "GT_CONF_PERCENTILE": 49.86,
                        "POS": 2289010,
                        "REF": "c",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "2289193c>t",
                    "nucleotide_index": 2289193,
                    "gene_name": "pncA",
                    "gene_position": 17,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 89.0,
                        "ALLELE_DP": [5.0, 84.0],
                        "FRS": 0.9438,
                        "COV_TOTAL": 89,
                        "COV": [5, 84],
                        "GT_CONF": 0.57,
                        "GT_CONF_PERCENTILE": 85.13,
                        "POS": 2289193,
                        "REF": "c",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "a1338c", "gene": "rrs", "gene_position": 1338},
                {"mutation": "a1401g", "gene": "rrs", "gene_position": 1401},
                {
                    "mutation": "T350S",
                    "gene": "rpoB",
                    "gene_position": 350,
                    "ref": "acc",
                    "alt": "tcc",
                },
                {
                    "mutation": "S450L",
                    "gene": "rpoB",
                    "gene_position": 450,
                    "ref": "tcg",
                    "alt": "ttg",
                },
                {
                    "mutation": "A90A",
                    "gene": "gyrA",
                    "gene_position": 90,
                    "ref": "gcg",
                    "alt": "gca",
                },
                {"mutation": "g270a", "gene": "gyrA", "gene_position": 270},
                {
                    "mutation": "S95T",
                    "gene": "gyrA",
                    "gene_position": 95,
                    "ref": "agc",
                    "alt": "acc",
                },
                {
                    "mutation": "G17S",
                    "gene": "pncA",
                    "gene_position": 17,
                    "ref": "ggc",
                    "alt": "agc",
                },
                {
                    "mutation": "G78C",
                    "gene": "pncA",
                    "gene_position": 78,
                    "ref": "ggc",
                    "alt": "tgc",
                },
                {
                    "gene": "fabG1",
                    "mutation": "L203X",
                    "gene_position": 203,
                    "ref": "ctg",
                    "alt": "xxt",
                },
            ],
            "effects": {
                "KAN": [
                    {
                        "gene": "rrs",
                        "mutation": "a1338c",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "rrs",
                        "mutation": "a1401g",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "RIF": [
                    {
                        "gene": "rpoB",
                        "mutation": "T350S",
                        "prediction": "U",
                        "evidence": {},
                    },
                    {
                        "gene": "rpoB",
                        "mutation": "S450L",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "MXF": [
                    {
                        "gene": "gyrA",
                        "mutation": "A90A",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "gyrA",
                        "mutation": "S95T",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "S"},
                ],
                "PZA": [
                    {
                        "gene": "pncA",
                        "mutation": "G17S",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "pncA",
                        "mutation": "G78C",
                        "prediction": "R",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
                "INH": [
                    {
                        "gene": "fabG1",
                        "mutation": "L203X",
                        "prediction": "F",
                        "evidence": {},
                    },
                    {"phenotype": "F"},
                ],
            },
            "antibiogram": {"KAN": "R", "RIF": "R", "MXF": "S", "PZA": "R", "INH": "F"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open("tests/outputs/15/resistance_prediction_report.json", "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))

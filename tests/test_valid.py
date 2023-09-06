'''Suite of unit tests based on the `gnomonicus` test suite
'''
import json

import gnomonicus
import pytest
#Helpful function import for testing nested JSON equality as it gives exact differences
from recursive_diff import recursive_eq

'''
Due to complications testing equalities of nested jsons of lists/dicts, there is a lot of 
code specificially dedicated to ensuring the lists are in the same order (they differ due to
a mixture of dictionary behaviour and different positions within files). However, we only care that
the contents of the JSON is there for these tests rather than caring about order.
'''

def prep_json(j: dict) -> dict:
    '''Prepare a JSON for comparison by removing fields which cannot be reproduced

    Args:
        j (dict): Initial JSON

    Returns:
        dict: JSON without fields such as time and file paths
    '''
    del j['meta']['time_taken_s']
    del j['meta']['UTC-datetime-completed']
    del j['meta']['catalogue_file']
    del j['meta']['reference_file']
    del j['meta']['vcf_file']
    return j

def variants_key(x):
    '''Used as the sorted(key=) function for reliably sorting the variants list

    Args:
        x (list): List of the ordered values as key, value pairs

    Returns:
        str: String of the `variant+gene_name`
    '''
    variant = ''
    gene = ''
    for i in x:
        if i[0] == 'variant':
            variant = i[1]
        elif i[0] == 'gene_name':
            gene = i[1]
    return variant+gene

def ordered(obj):
    '''Recursively sort a JSON for equality checking. Based on https://stackoverflow.com/questions/25851183/how-to-compare-two-json-objects-with-the-same-elements-in-a-different-order-equa

    Args:
        obj (object): Any JSON element. Probably one of dict, list, tuple, str, int, None

    Returns:
        object: Sorted JSON
    '''
    if isinstance(obj, dict):
        if 'variants' in obj.keys():
            #We have the 'data' field which needs a little extra nudge
            if 'effects' in obj.keys():
                #Case when we have effects populated
                return [
                    ('antibiogram', ordered(obj['antibiogram'])),
                    ('effects', ordered(obj['effects'])),
                    ('mutations', ordered(obj['mutations'])),
                    ('variants', sorted([ordered(x) for x in obj['variants']], key=variants_key))  
                ]
            elif 'mutations' in obj.keys():
                #Case for if we have mutations and variants but no effects
                return [
                    ('mutations', ordered(obj['mutations'])),
                    ('variants', sorted([ordered(x) for x in obj['variants']], key=variants_key))
                ]                
            else:
                #Case where no effects/mutations have been populated
                return [
                    ('variants', sorted([ordered(x) for x in obj['variants']], key=variants_key))
                ]
        else:
            return sorted((k, ordered(obj[k])) for k in sorted(list(obj.keys())))

    if isinstance(obj, list) or isinstance(obj, tuple):
        return sorted(ordered(x) for x in obj)
    
    #Nones cause issues with ordering as there is no < operator.
    #Convert to string to avoid this
    if type(obj) == type(None):
        return str(obj)
    
    #Because nan types are helpful, `float('nan') == float('nan') -> False`
    #So check based on str value rather than direct equality
    if str(obj) == 'nan':
        #Conversion to None for reproducability
        return str(None)
    
    if isinstance(obj, int):
        #Ints are still ordered (just not numerically) if sorted by str value, so convert to str
        #to allow sorting lists of None/int
        return str(obj)
    if isinstance(obj, float):
        #Similarly convert float, but check if they are x.0 to avoid str comparison issues
        if int(obj) == obj:
            return str(int(obj))
        else:
            return str(obj)
    else:
        return obj

def test_1():
    '''Input:
            NC_045512.2-S_E484K-minos.vcf
        Expect output:
            variants:    23012g>a
            mutations:   S@E484K
            predictions: {'AAA': 'R', 'BBB': 'S'}
    '''
    vcfStem = "NC_045512.2-S_E484K-minos"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '23012g>a',
                    'nucleotide_index': 23012,
                    'gene_name': 'S',
                    'gene_position': 484,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        'GT': [1, 1], 'DP': 44, 'DPF': 0.991, 'COV': [0, 44], 
                        'FRS': 1.0, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 
                        'REF': 'g', 'ALTS': ['a'], 'POS': 23012
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'E484K',
                    'gene': 'S',
                    'gene_position':484,
                    "ref": "gaa",
                    "alt": "aaa"
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'E484K',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/1/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))

def test_2():
    '''Input:
            NC_045512.2-S_E484K-samtools.vcf
        Expect output:
            variants:    23012g>a
            mutations:   S@E484K
            predictions: {'AAA': 'R', 'BBB': 'S'}
    '''
    vcfStem = "NC_045512.2-S_E484K-samtools"
    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '23012g>a',
                    'nucleotide_index': 23012,
                    'gene_name': 'S',
                    'gene_position': 484,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "PL": [
                            255,
                            33,
                            0
                        ],
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": [
                            "a"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'E484K',
                    'gene': 'S',
                    'gene_position':484,
                    "ref": "gaa",
                    "alt": "aaa"
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'E484K',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/2/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))

def test_3():
    '''Input:
            NC_045512.2-S_F2F-minos.vcf
        Expect output:
            variants:    21568t>c
            mutations:   S@F2F
            predictions: {'AAA': 'S', 'BBB': 'S'}
    '''
    vcfStem = "NC_045512.2-S_F2F-minos"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '21568t>c',
                    'nucleotide_index': 21568,
                    'gene_name': 'S',
                    'gene_position': 2,
                    'codon_idx': 2,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21568,
                        "REF": "t",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'F2F',
                    'gene': 'S',
                    'gene_position': 2,
                    'ref': 'ttt',
                    'alt': 'ttc'
                },
                {
                    'mutation': 't6c',
                    'gene': 'S',
                    'gene_position': 6
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'F2F',
                        'prediction': 'S',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'S'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'S',
                'BBB': 'S'
            }
        }
    }


    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/3/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_4():
    '''Input:
            NC_045512.2-S_F2L-minos.vcf
        Expect output:
            variants:    21566t>c
            mutations:   S@F2L
            predictions: {'AAA': 'U', 'BBB': 'S'}
    '''
    #Setup
    vcfStem = "NC_045512.2-S_F2L-minos"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '21566t>c',
                    'nucleotide_index': 21566,
                    'gene_name': 'S',
                    'gene_position': 2,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21566,
                        "REF": "t",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'F2L',
                    'gene': 'S',
                    'gene_position': 2,
                    'ref': 'ttt',
                    'alt': 'ctt'
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'F2L',
                        'prediction': 'U',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'U'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'U',
                'BBB': 'S'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/4/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))
   


def test_5():
    '''Input:
            NC_045512.2-S_200_indel-minos.vcf
        Expect output:
            variants:    21762_ins_c
            mutations:   S@200_ins_c
            predictions: {'AAA': 'R', 'BBB': 'S'}
    '''
    vcfStem = "NC_045512.2-S_200_indel-minos"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '21762_ins_c',
                    'nucleotide_index': 21762,
                    'gene_name': 'S',
                    'gene_position': 200,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21762,
                        "REF": "c",
                        "ALTS": [
                            "cc"
                        ]
                    },
                    'vcf_idx': 1
                },
            ],
            'mutations': [
                {
                    'mutation': '200_ins_c',
                    'gene': 'S',
                    'gene_position':200
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': '200_ins_c',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/5/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_6():
    '''Input:
            NC_045512.2-double-minos.vcf
        Expect output:
            variants:    27758g>c
            mutations:   ORF7a!122S, ORF7b@M1I
            predictions: {'AAA': 'R', 'BBB': 'R'}
    '''
    vcfStem = "NC_045512.2-double-minos"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '27758g>c',
                    'nucleotide_index': 27758,
                    'gene_name': 'ORF7a',
                    'gene_position': 122,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '27758g>c',
                    'nucleotide_index': 27758,
                    'gene_name': 'ORF7b',
                    'gene_position': 1,
                    'codon_idx': 2,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': '!122S',
                    'gene': 'ORF7a',
                    'gene_position': 122,
                    'ref': 'tga',
                    'alt': 'tca'
                },
                {
                    'mutation': 'M1I',
                    'gene': 'ORF7b',
                    'gene_position': 1,
                    'ref': 'atg',
                    'alt': 'atc'
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'ORF7a',
                        'mutation': '!122S',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
                'BBB': [
                    {
                        'gene': 'ORF7b',
                        'mutation': 'M1I',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'R'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/6/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_7():
    '''Input:
            NC_045512.2-S_E484K&1450_ins_a-minos.vcf
        Expect output:
            variants:    23012g>a, 23012_ins_a
            mutations:   S@E484K, S@1450_ins_a, S@1450_ins_a&S@E484K
            predictions: {'AAA': 'R', 'BBB': 'R'}
    '''
    vcfStem = "NC_045512.2-S_E484K&1450_ins_a-minos"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '23012_ins_a',
                    'nucleotide_index': 23012,
                    'gene_name': 'S',
                    'gene_position': 1450,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": [
                            "aa"
                        ]
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '23012g>a',
                    'nucleotide_index': 23012,
                    'gene_name': 'S',
                    'gene_position': 484,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": [
                            "aa"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': '1450_ins_a',
                    'gene': 'S',
                    'gene_position':1450
                },
                {
                    'mutation': 'E484K',
                    'gene': 'S',
                    'gene_position':484,
                    'ref': "gaa",
                    "alt": "aaa"
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'E484K',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'gene': 'S',
                        'mutation': '1450_ins_a',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
                'BBB': [
                    {
                        'gene': None,
                        'mutation': 'S@1450_ins_a&S@E484K',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'R'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/7/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))

def test_8():
    '''Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:0.045, 25283_del_g:0.045, 25252_ins_cc:0.045
        mutations:   !1274Q:0.045, 3721_del_g:0.045, 3690_ins_cc:0.045
        predictions: {'AAA': 'R'}
    '''
    vcfStem = "NC_045512.2-minors"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '25382t>c:0.045',
                    'nucleotide_index': 25382,
                    'gene_name': 'S',
                    'gene_position': 1274,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25382,
                        "REF": "t",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '21558g>a:0.045',
                    'nucleotide_index': 21558,
                    'gene_name': 'S',
                    'gene_position': -5,
                    'codon_idx': None,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21558,
                        "REF": "g",
                        "ALTS": [
                            "a"
                        ]
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '25252_ins_cc:0.045',
                    'nucleotide_index': 25252,
                    'gene_name': 'S',
                    'gene_position': 3690,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25252,
                        "REF": "g",
                        "ALTS": [
                            "gcc"
                        ]                         
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '25283_del_t:0.045',
                    'nucleotide_index': 25283,
                    'gene_name': 'S',
                    'gene_position': 3721,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25282,
                        "REF": "tt",
                        "ALTS": [
                            "t"
                        ]                       
                    },
                    'vcf_idx': 1
                },                
            ],
            'mutations': [
                {
                    'mutation': '!1274Q:0.045',
                    'gene': 'S',
                    'gene_position':1274,
                    'ref': 'taa',
                    'alt': 'zzz'
                },
                {
                    'mutation': 'g-5a:0.045',
                    'gene': 'S',
                    'gene_position':-5
                },
                {
                    'mutation': '3721_del_t:0.045',
                    'gene': 'S',
                    'gene_position':3721
                },
                {
                    'mutation': '3690_ins_cc:0.045',
                    'gene': 'S',
                    'gene_position':3690
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': '!1274Q:0.045',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'gene': 'S',
                        'mutation': 'g-5a:0.045',
                        'prediction': 'U',
                        'evidence': {}
                    },
                    {
                        'gene': 'S',
                        'mutation': '3721_del_t:0.045',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'gene': 'S',
                        'mutation': '3690_ins_cc:0.045',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/8/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_9():
    '''Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:2, 25283_del_g:2, 25252_ins_cc:2
        mutations:   !1274Q:2, 3721_del_g:2, 3690_ins_cc:2
        predictions: {'AAA': 'R'}
    '''
    vcfStem = "NC_045512.2-minors"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '25382t>c:2',
                    'nucleotide_index': 25382,
                    'gene_name': 'S',
                    'gene_position': 1274,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25382,
                        "REF": "t",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '21558g>a:2',
                    'nucleotide_index': 21558,
                    'gene_name': 'S',
                    'gene_position': -5,
                    'codon_idx': None,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21558,
                        "REF": "g",
                        "ALTS": [
                            "a"
                        ]
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '25252_ins_cc:2',
                    'nucleotide_index': 25252,
                    'gene_name': 'S',
                    'gene_position': 3690,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25252,
                        "REF": "g",
                        "ALTS": [
                            "gcc"
                        ]                         
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '25283_del_t:2',
                    'nucleotide_index': 25283,
                    'gene_name': 'S',
                    'gene_position': 3721,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            0,
                            0
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            42,
                            2
                        ],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25282,
                        "REF": "tt",
                        "ALTS": [
                            "t"
                        ]                       
                    },
                    'vcf_idx': 1
                },                
            ],
            'mutations': [
                {
                    'mutation': '!1274Q:2',
                    'gene': 'S',
                    'gene_position':1274,
                    'ref': 'taa',
                    'alt': 'zzz'
                },
                {
                    'mutation': 'g-5a:2',
                    'gene': 'S',
                    'gene_position':-5
                },
                {
                    'mutation': '3721_del_t:2',
                    'gene': 'S',
                    'gene_position':3721
                },
                {
                    'mutation': '3690_ins_cc:2',
                    'gene': 'S',
                    'gene_position':3690
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': '!1274Q:2',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'gene': 'S',
                        'mutation': 'g-5a:2',
                        'prediction': 'U',
                        'evidence': {}
                    },
                    {
                        'gene': 'S',
                        'mutation': '3721_del_t:2',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'gene': 'S',
                        'mutation': '3690_ins_cc:2',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/9/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_10():
    '''Testing a catalogue and sample which have large deletions
    Input:
        TEST-DNA-large-del.vcf
    Expect output:
        variants:    3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc
        mutations:   A@-1_del_aaaaaaaaccccccccccgggggggggg, A@del_0.93, B@del_1.0, C@4_del_ggg
        predictions: {'AAA': 'R'}
    '''
    vcfStem = "TEST-DNA-large-del"

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "TEST_DNA",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test_dna",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc',
                    'nucleotide_index': 3,
                    'gene_name': 'A',
                    'gene_position': -1,
                    'codon_idx': None,
                        "vcf_evidence": {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 2,
                        "COV": [
                            1,
                            1
                        ],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": [
                            "a"
                        ]
                        },
                    'vcf_idx': 1
                },
                {
                    'variant': '3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc',
                    'nucleotide_index': 3,
                    'gene_name': 'B',
                    'gene_position': 1,
                    'codon_idx': 0,
                        "vcf_evidence": {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 2,
                        "COV": [
                            1,
                            1
                        ],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": [
                            "a"
                        ]
                        },
                    'vcf_idx': 1
                },
                {
                    'variant': '3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc',
                    'nucleotide_index': 3,
                    'gene_name': 'C',
                    'gene_position': 4,
                    'codon_idx': 2,
                        "vcf_evidence": {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 2,
                        "COV": [
                            1,
                            1
                        ],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": [
                            "a"
                        ]
                        },
                    'vcf_idx': 1
                },
            ],
            'mutations': [
                {
                    'mutation': '-1_del_aaaaaaaaccccccccccgggggggggg',
                    'gene': 'A',
                    'gene_position':-1,
                },
                {
                    'mutation': 'del_0.93',
                    'gene': 'A',
                    'gene_position': None
                },
                {
                    'mutation': 'del_1.0',
                    'gene': 'B',
                    'gene_position': None
                },
                {
                    'mutation': '4_del_ggg',
                    'gene': 'C',
                    'gene_position': 4
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'A',
                        'mutation': '-1_del_aaaaaaaaccccccccccgggggggggg',
                        'prediction': 'U',
                        'evidence': {}
                    },
                    {
                        'gene': 'A',
                        'mutation': 'del_0.93',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'gene': 'B',
                        'mutation': 'del_1.0',
                        'prediction': 'U',
                        'evidence': {}
                    },
                    {
                        'gene': 'C',
                        'mutation': '4_del_ggg',
                        'prediction': 'U',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open("tests/outputs/10/resistance_prediction_report.json", 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))

#TODO: Add sythetic samples for MDR/XDR etc (and maybe a real sample?)

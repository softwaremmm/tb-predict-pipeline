import json
import os

import pandas as pd
import pytest
from recursive_diff import recursive_eq


def concatFields(d: dict) -> str:
    '''Concat the value of a dictionary in the order of the sorted keys

    Args:
        d (dict): Dictionary input

    Returns:
        str: String of values
    '''
    return ''.join([str(d[key]) for key in sorted(list(d.keys()))])

def sortValues(json: dict) -> dict:
    '''Sort the values within the VARIANTS, MUTATIONS and EFFECTS lists in a JSON.
    THis allows to test for contents equality as order is not particularly important here

    Args:
        json (dict): JSON in

    Returns:
        dict: JSON with VARIANTS, MUTATIONS and EFFECTS lists in reproducable orders for equality checks
    '''
    variants = json['data']['VARIANTS']
    mutations = json['data'].get('MUTATIONS', None)
    effects = json['data'].get('EFFECTS', None)
    
    json['data']['VARIANTS'] = sorted(variants, key=concatFields)
    if mutations is not None:
        json['data']['MUTATIONS'] = sorted(mutations, key=concatFields)
    if effects is not None:
        for drug in effects.keys():
            json['data']['EFFECTS'][drug] = sorted(effects[drug], key=concatFields)
    
    return json

def test_1():
    '''A test case from the `gnomonicus` unit tests to ensure consistency.
    Input:
            NC_045512.2-double-minos.vcf
        Expect output:
            variants:    27758g>c
            mutations:   ORF7a!122S, ORF7b@M1I
            predictions: {'AAA': 'R', 'BBB': 'R'}
    '''

    path = "tests/outputs/1/NC_045512/"
    vcfStem = "NC_045512"

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '27758g>c'


    assert 'ORF7a' in mutations['GENE'].to_list()
    assert 'ORF7b' in mutations['GENE'].to_list()

    assert mutations['MUTATION'][mutations['GENE'].to_list().index('ORF7a')] == '!122S'
    assert mutations['MUTATION'][mutations['GENE'].to_list().index('ORF7b')] == 'M1I'

    assert 'AAA' in effects['DRUG'].to_list()
    assert 'BBB' in effects['DRUG'].to_list()
    
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('BBB')] == 'R'

    expectedJSON = {
        'meta': {
            'version': '1.1.3',
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ], 
                        "BBB": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ], 
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '27758g>c',
                    'NUCLEOTIDE_INDEX': 27758
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': '!122S',
                    'GENE': 'ORF7a',
                    'GENE_POSITION': 122
                },
                {
                    'MUTATION': 'M1I',
                    'GENE': 'ORF7b',
                    'GENE_POSITION': 1
                },
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'ORF7a',
                        'MUTATION': '!122S',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
                'BBB': [
                    {
                        'GENE': 'ORF7b',
                        'MUTATION': 'M1I',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)

def test_2():
    '''A test case which has a single variant, no mutations and no effects
    Input:
            NC_045512.2-just-variant.vcf
        Expect output:
            variants:    4a>g
    '''

    path = "tests/outputs/2/NC_045512/"
    vcfStem = "NC_045512"

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    with pytest.raises(Exception):
        _ = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    with pytest.raises(Exception):
        _ = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '4a>g'

    expectedJSON = {
        'meta': {
            'version': '1.1.3',
            'guid': vcfStem,
            'fields': {
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '4a>g',
                    'NUCLEOTIDE_INDEX': 4
                }
            ],
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)



def test_3():
    '''Testing a case which generates variants and mutations but not effects
    Input:
            NC_045512.2-double-minos.vcf
        Expect output:
            variants:    28280g>t
            mutations:   N@D3Y
    '''

    path = "tests/outputs/3/NC_045512/"
    vcfStem = "NC_045512"

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    with pytest.raises(Exception):
        _ = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '28280g>t'


    assert 'N' in mutations['GENE'].to_list()

    assert mutations['MUTATION'][mutations['GENE'].to_list().index('N')] == 'D3Y'

    

    expectedJSON = {
        'meta': {
            'version': '1.1.3',
            'guid': vcfStem,
            'fields': {
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '28280g>t',
                    'NUCLEOTIDE_INDEX': 28280
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': 'D3Y',
                    'GENE': 'N',
                    'GENE_POSITION': 3
                },
            ],
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)

def test_4():
    '''Test to ensure that using docker does not change outputs
    '''
    #Docker generated files
    path = "tests/outputs/1/NC_045512/"
    vcfStem = "NC_045512"

    #Check for expected values within csvs
    variants1 = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations1 = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects1 = pd.read_csv(path + f"{vcfStem}.effects.csv")

    JSON1 = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del JSON1['meta']['UTC-datetime-run']

    #Bare metal generated files
    path = "tests/outputs/4/NC_045512/"

    #Check for expected values within csvs
    variants2 = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations2 = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects2 = pd.read_csv(path + f"{vcfStem}.effects.csv")

    JSON2 = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del JSON2['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    print(json.dumps(JSON1, indent=2, sort_keys=True))
    print()
    print(json.dumps(JSON2, indent=2, sort_keys=True))
    recursive_eq(JSON1, JSON2)
    recursive_eq(variants1, variants2)
    recursive_eq(mutations1, mutations2)
    recursive_eq(effects1, effects2)
"""Comapre two gnomonicus JSONs for results"""
import json
import sys
import argparse


def compare_meta(j1: dict, j2: dict) -> int:
    """Compare the meta data of the two JSONs"""
    diff = 0
    m1 = j1["meta"]
    m2 = j2["meta"]

    relevant_keys = [
        "status",
        "workflow_version",
        "reference",
        "catalogue_type",
        "catalogue_name",
        "catalogue_version",
    ]
    for key in relevant_keys:
        if m1[key] != m2[key]:
            diff += 1
            print(f"\t{key}: {m1[key]} != {m2[key]}")
    return diff


def compare_antibiogram(j1: dict, j2: dict) -> int:
    """Compare the antibiogram data of the two JSONs"""
    diff = 0
    a1 = j1["data"]["antibiogram"]
    a2 = j2["data"]["antibiogram"]
    drugs1 = set(a1.keys())
    drugs2 = set(a2.keys())
    extra1 = drugs1 - drugs2
    extra2 = drugs2 - drugs1
    if extra1:
        diff += 1
        print(f"\tLHS has extra antibiogram drugs: {extra1}")
    if extra2:
        diff += 1
        print(f"\tRHS has extra antibiogram drugs: {extra2}")
    for drug in sorted(drugs1 & drugs2):
        if a1[drug] != a2[drug]:
            diff += 1
            print(f"\t{drug}: {a1[drug]} != {a2[drug]}")

    return diff


def compare_variants(j1: dict, j2: dict) -> int:
    """Compare the variants of the two JSONs"""
    diff = 0
    vs1: dict[str, str] = {}
    vs2: dict[str, str] = {}

    # Need to turn variant lists into easier to compare format
    for v in j1["data"]["variants"]:
        key = f"{v['gene_name']}{v['gene_position']}"
        value = v["variant"]
        vs1[key] = value

    for v in j2["data"]["variants"]:
        key = f"{v['gene_name']}{v['gene_position']}"
        value = v["variant"]
        vs2[key] = value

    all_keys = set(vs1.keys()) | set(vs2.keys())
    for key in sorted(list(all_keys)):
        if key not in vs2:
            # Must be in vs1 then
            diff += 1
            print(f"\t{key}: {vs1[key]} != missing")
            continue
        if key not in vs1:
            diff += 1
            print(f"\t{key}: missing != {vs2[key]}")
            continue

        if vs1[key] != vs2[key]:
            diff += 1
            print(f"\t{key}: {vs1[key]} != {vs2[key]}")

    return diff


def compare_mutations(j1: dict, j2: dict) -> int:
    """Compare the mutations of the two JSONs"""
    diff = 0
    muts1: dict[str, str] = {}
    muts2: dict[str, str] = {}

    # Need to turn variant lists into easier to compare format
    for mut in j1["data"]["mutations"]:
        key = f"{mut['gene']}{mut['gene_position']}"
        value = mut["mutation"]
        muts1[key] = value

    for mut in j2["data"]["mutations"]:
        key = f"{mut['gene']}{mut['gene_position']}"
        value = mut["mutation"]
        muts2[key] = value

    all_keys = set(muts1.keys()) | set(muts2.keys())
    for key in sorted(list(all_keys)):
        if key not in muts2:
            # Must be in vs1 then
            diff += 1
            print(f"\t{key}: {muts1[key]} != missing")
            continue
        if key not in muts1:
            diff += 1
            print(f"\t{key}: missing != {muts2[key]}")
            continue

        if muts1[key] != muts2[key]:
            diff += 1
            print(f"\t{key}: {muts1[key]} != {muts2[key]}")

    return diff


def get_effect_causes(drug_effect_list: list[dict]) -> tuple[str, set[str]]:
    """Get the phenotype and causes from a drug effect list"""
    phenotype = ""
    causes = []
    for cause in drug_effect_list:
        if "phenotype" in cause:
            phenotype = cause["phenotype"]
            continue
        cause_str = f"{cause['gene']}:{cause['mutation']} ({cause['prediction']})"
        causes.append(cause_str)

    if phenotype == "":
        raise ValueError(
            f"Effects should always include a phenotype prediction. {drug_effect_list}"
        )
    return (phenotype, set(causes))


def compare_effects(j1: dict, j2: dict) -> int:
    """Compare the effects of the two JSONs"""
    diff = 0

    effects1 = j1["data"]["effects"]
    effects2 = j2["data"]["effects"]
    all_drugs = set(effects1.keys()) | set(effects2.keys())

    for drug in sorted(all_drugs):
        if drug not in effects2:
            # must be in drugs 1
            diff += 1
            phenotype1, causes1 = get_effect_causes(effects1[drug])
            print(f"effects-{drug}: {phenotype1} != No effect (S)")
            print("\tCauses for LHS effect:")
            for cause in sorted(causes1):
                print(f"\t\t{cause}")
            continue
        if drug not in effects1:
            diff += 1
            phenotype2, causes2 = get_effect_causes(effects2[drug])
            print(f"effects-{drug}: No effect (S) != {phenotype2}")
            print("\tCauses for RHS effect:")
            for cause in sorted(causes2):
                print(f"\t\t{cause}")
            continue

        phenotype1, causes1 = get_effect_causes(effects1[drug])
        phenotype2, causes2 = get_effect_causes(effects2[drug])
        drug_diff = False
        if phenotype1 != phenotype2:
            diff += 1
            drug_diff = True
            print(f"effects-{drug}:")
            print(f"\tphenotype: {phenotype1} != {phenotype2}")

        for cause in sorted(causes1 ^ causes2):
            if not drug_diff:
                drug_diff = True
                print(f"effects-{drug}:")
                print(f"\tphenotype: {phenotype1} == {phenotype2}")

            if cause in causes1:
                diff += 1
                print(f"\t{cause} only in LHS")
            else:
                diff += 1
                print(f"\t{cause} only in RHS")

    return diff


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


def compare(path1: str, path2: str) -> None:
    """Load the 2 JSONs and compare for equality

    Args:
        path1 (str): Path to JSON 1
        path2 (str): Path to JSON 2
    """
    json1 = prep_json(json.load(open(path1, "r", encoding="utf-8")))
    json2 = prep_json(json.load(open(path2, "r", encoding="utf-8")))

    diff = 0
    print("---- Meta ----")
    diff += compare_meta(json1, json2)
    print("---- Variants ----")
    diff += compare_variants(json1, json2)
    print("---- Mutations ----")
    diff += compare_mutations(json1, json2)
    print("---- Effects ----")
    diff += compare_effects(json1, json2)
    print("---- Antibiogram ----")
    diff += compare_antibiogram(json1, json2)

    print(f"Found {diff} differences")
    if diff > 0:
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare two gnomonicus JSONs for results"
    )
    parser.add_argument("json1", help="Path to JSON 1")
    parser.add_argument("json2", help="Path to JSON 2")
    args = parser.parse_args()
    compare(args.json1, args.json2)

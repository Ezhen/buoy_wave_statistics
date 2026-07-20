"""
Builds stage_registry.json by scanning every numbered pipeline stage
script and parsing ITS OWN --min-years argparse default directly from
source, via AST - not by hand-maintaining a mirrored dict.

This closes a real maintenance gap: build_buoy_registry.py's
STAGE_MIN_YEARS dict hardcodes 10.0/3.0/10.0 for Stages 14/15/25 and
says so explicitly in its own docstring - "if any of these stage
defaults change, update STAGE_MIN_YEARS below to match... this file
does not read the other scripts' source at runtime, so it can silently
drift out of sync." This script IS that runtime read, so the value is
always current with whatever the stage script itself actually does,
with no separate value to remember to update.

Uses AST parsing rather than importing/executing each stage script
(which would require its dependencies installed and could have side
effects) or regex (fragile against formatting variation) - this reads
the parser.add_argument(...) call structure directly.

Also reuses run_all_buoys.py's ALREADY-EXISTING requirements schema
(variables_any/variables_all/min_record_years, consumed by its
stage_eligible() function) rather than inventing a new one - Stage 09's
variables_any requirement is carried over from there directly, since
that's the one variable-requirement case already confirmed correct and
in active use; variable requirements for stages 14+ are NOT
auto-derivable from source the way --min-years is (they're baked into
each script's internal logic, not a CLI default), so this section is
manually curated and explicitly marked as such, same pattern as
build_buoy_registry.py's KNOWN_QUIRKS.

Usage:
    python tools/build_stage_registry.py
"""

import ast
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Manually curated - NOT auto-derived, unlike min_record_years below.
# Carried over from run_all_buoys.py's existing STAGES list, the one
# already-confirmed-correct case, plus scripts checked directly this
# session where the requirement is unambiguous from source (e.g. Stage
# 16's VMDR-only directional-alignment check is conditional, not a hard
# requirement for the whole stage, so deliberately left out here rather
# than guessed at).
KNOWN_VARIABLE_REQUIREMENTS = {
    "09_cross_variable_analysis.py": {"variables_any": ["VTPK", "VMDR"]},
}


def extract_min_years_default(script_path: Path):
    """AST-parse a stage script's own --min-years argparse default.
    Returns None if the script has no --min-years argument at all -
    that's a real, meaningful result (most stages don't gate on record
    length), not a parsing failure."""
    try:
        tree = ast.parse(script_path.read_text())
    except SyntaxError:
        return None, "syntax error - could not parse"

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        is_add_argument = isinstance(func, ast.Attribute) and func.attr == "add_argument"
        if not is_add_argument:
            continue
        if not node.args or not isinstance(node.args[0], ast.Constant):
            continue
        if node.args[0].value != "--min-years":
            continue
        for kw in node.keywords:
            if kw.arg == "default" and isinstance(kw.value, ast.Constant):
                return kw.value.value, None
        return None, "--min-years found but default isn't a simple literal - check manually"

    return None, None  # no --min-years argument in this script at all


def main():
    # Require a 2-digit numeric prefix (02, 03b, 11b, 12b, ... 26) - a
    # single leading digit like "0_download_..." is the one-time data-
    # fetch utility, not a per-buoy analysis stage with the same
    # --buoy/--var invocation contract as everything else. Including it
    # here would be a real wrong inclusion once this registry drives a
    # generic runner that iterates "every stage, for every buoy."
    stage_scripts = sorted(
        p for p in REPO_ROOT.glob("*.py")
        if p.stem[:2].isdigit() and p.name != "01_load_clean.py"
    )

    registry = {}
    parse_issues = []
    for script in stage_scripts:
        min_years, issue = extract_min_years_default(script)
        if issue:
            parse_issues.append(f"{script.name}: {issue}")

        requirements = {}
        if min_years is not None:
            requirements["min_record_years"] = min_years
        if script.name in KNOWN_VARIABLE_REQUIREMENTS:
            requirements.update(KNOWN_VARIABLE_REQUIREMENTS[script.name])

        registry[script.name] = {
            "min_record_years": min_years,
            "requirements": requirements,
            "has_variable_requirement_curated": script.name in KNOWN_VARIABLE_REQUIREMENTS,
        }

    out_path = REPO_ROOT / "stage_registry.json"
    with open(out_path, "w") as f:
        json.dump(registry, f, indent=2)

    print(f"Wrote {out_path}: {len(registry)} stage script(s) scanned.")
    n_with_gate = sum(1 for v in registry.values() if v["min_record_years"] is not None)
    print(f"  {n_with_gate} script(s) have a --min-years gate:")
    for name, entry in sorted(registry.items()):
        if entry["min_record_years"] is not None:
            print(f"    {name}: {entry['min_record_years']}yr")
    print(f"  {len(KNOWN_VARIABLE_REQUIREMENTS)} script(s) have a manually-curated variable "
          f"requirement (see KNOWN_VARIABLE_REQUIREMENTS in this file - NOT auto-derived, "
          f"unlike min_record_years above).")
    if parse_issues:
        print(f"\n  Parse issues found (review manually):")
        for issue in parse_issues:
            print(f"    {issue}")


if __name__ == "__main__":
    main()

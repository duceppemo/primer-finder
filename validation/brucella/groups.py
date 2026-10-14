"""Which species make up the inclusion and exclusion group of each target.

Four scripts of this harness need the same answer and they are not a package, so the definitions live
here and each of them imports it. A target named here is a question about the panel, not a setting of
primer-finder:

- `genus`: the eight classical species the paper's in-silico panel holds, against the non-classical
  relatives. This is the paper's own question and is never changed, so its record stays reproducible;
- `genus_wide`: the same plus the other four species of the classical clade -- B. pinnipedialis,
  B. amazoniensis, B. vulpis and B. inopinata -- which the paper's panel has none of. What a genus assay
  would have to cover to span the whole clade;
- `abortus`, `melitensis`, `suis`: one species against the other seven of the paper's eight, and the
  relatives;
- `suis_lato`: **B. suis together with B. canis**, against the other six of the eight. B. canis is a
  lineage *inside* B. suis phylogenetically, so no stretch is present in every B. suis and absent from
  every B. canis: the `suis` target reports nothing at any setting, and this is the same question asked
  in a way the phylogeny allows. An assay from it detects both, which is what "B. suis sensu lato" means.

The four `classical_extra` species are in neither group of a species target, as they are in neither group
of `genus`: the paper's panel has none of them, and putting them in one would change what the groups mean.
`Brucella sp.` is in neither group of any target, because roughly half of those assemblies are former
Ochrobactrum.
"""

from __future__ import annotations

# The eight classical Brucella species the paper's in-silico panel holds (its Table S2).
CLASSICAL = ("abortus", "melitensis", "suis", "canis", "ceti", "ovis", "neotomae", "microti")

# The other species of the classical clade, which the paper's panel does not hold.
CLASSICAL_EXTRA = ("pinnipedialis", "amazoniensis", "vulpis", "inopinata")

# Targets that are not a single species: what their inclusion group is made of.
COMPOSITE = {
    "genus": CLASSICAL,
    "genus_wide": CLASSICAL + CLASSICAL_EXTRA,
    "suis_lato": ("suis", "canis"),
}

TARGETS = ("genus", "genus_wide", "abortus", "melitensis", "suis", "suis_lato")


def inclusion_species(target: str) -> tuple[str, ...]:
    """The species the inclusion group of a target is made of."""
    return COMPOSITE.get(target, (target,))


def exclusion_species(target: str) -> tuple[str, ...]:
    """The classical species its exclusion group is made of, besides the non-classical relatives.

    A genus target has none: everything classical is in its inclusion group. Any other target is set
    against the species of the paper's eight that are not in its own inclusion group.
    """
    if target.startswith("genus"):
        return ()
    mine = set(inclusion_species(target))
    return tuple(name for name in CLASSICAL if name not in mine)

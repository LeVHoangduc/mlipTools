def get_E0(db, key, config_type="IsolatedAtom"):
    """Extract isolated-atom reference energies."""

    E0 = {}

    for atoms in db:
        if (
            atoms.info.get("config_type") == config_type
            and len(atoms) == 1
        ):
            element = atoms.get_chemical_symbols()[0]
            E0[element] = float(atoms.info[key])

    return E0


def get_atomization_energies(
    atoms,
    E0,
    energy_key,
    peratom=False,
):
    """Calculate atomization energy."""

    E_total = float(atoms.info[energy_key])
    symbols = atoms.get_chemical_symbols()

    missing = sorted(set(symbols) - set(E0))

    if missing:
        raise KeyError(
            "Missing isolated-atom energies for: "
            + ", ".join(missing)
        )

    E_isolated = sum(E0[element] for element in symbols)
    E_atomization = E_total - E_isolated

    if peratom:
        E_atomization /= len(atoms)

    return E_atomization
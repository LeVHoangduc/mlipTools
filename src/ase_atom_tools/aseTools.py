import numpy as np

def get_E0s(db, key, config_type="IsolatedAtom"):
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


def get_atomization_energy(
    atoms,
    E0,
    energy_key,
    peratom=False,
):
    """Get atomization energy for each configuration"""

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


def get_total_energy(db, energy_key, peratom=False):
    """Get total energy or total enery per atom for a database, shape (nconfig,)"""

    energies = []
    for atoms in db:
        E = float(atoms.info[energy_key])
        if peratom:
            E /= len(atoms)
        energies.append(E)

    return np.array(energies)

def get_forces(
        db,
        forces_key,
):
    """Get forces for a database: list of nconfig arrays, each of shape (natoms_i, 3)"""

    return [atoms.arrays[forces_key] for atoms in db]

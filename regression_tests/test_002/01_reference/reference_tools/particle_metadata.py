"""Write particle metadata for Test 002 reference files."""

import h5py

from regression_test_tools import ConfigError


# Physical isotope masses match the values used by FIDASIM.
SPECIES_MASS_AMU = {
    "h": 1.007276466879,
    "d": 2.013553212745,
    "t": 3.01550071632,
}


def write_particle_metadata(h5file, particle):
    """Write the scalar particle-metadata datasets to a reference HDF5 file."""
    species = str(particle["species"]).strip().lower()
    if species not in SPECIES_MASS_AMU:
        raise ConfigError(f"Unsupported particle species: {species}")

    string_type = h5py.string_dtype(encoding="utf-8")
    species_dataset = h5file.create_dataset(
        "species",
        data=species,
        dtype=string_type,
    )
    atomic_number = h5file.create_dataset(
        "atomic_number",
        data=int(particle["atomic_number"]),
    )
    mass_number = h5file.create_dataset(
        "mass_number",
        data=int(particle["mass_number"]),
    )
    charge_state = h5file.create_dataset(
        "charge_state",
        data=int(particle["charge_state"]),
    )
    isotope_mass = h5file.create_dataset(
        "A",
        data=SPECIES_MASS_AMU[species],
    )

    species_dataset.attrs["description"] = "Canonical fast-ion species identifier"
    atomic_number.attrs["description"] = "Number of protons in the ion nucleus"
    mass_number.attrs["description"] = (
        "Integer isotope mass number: protons plus neutrons"
    )
    charge_state.attrs["description"] = (
        "Ion charge state in elementary-charge units"
    )
    isotope_mass.attrs["description"] = "Physical isotope mass"

    atomic_number.attrs["units"] = "dimensionless"
    mass_number.attrs["units"] = "dimensionless"
    charge_state.attrs["units"] = "elementary charge"
    isotope_mass.attrs["units"] = "amu"

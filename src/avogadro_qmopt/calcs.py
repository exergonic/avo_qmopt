import logging
import re

import psi4

from . import CALCS_DIR

BOHR_TO_ANGSTROM = 0.52917721090380

logger = logging.getLogger(__name__)


_ELEM_SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
]


def _sanitize_name(name: str) -> str:
    return re.sub(r"[^\w-]", "_", name.strip()) or "molecule"


def _option(options: dict, key: str, default):
    return options.get(key, default)


def _get_elem_symbol(z: int) -> str:
    return _ELEM_SYMBOLS[z] if 0 <= z < len(_ELEM_SYMBOLS) else "X"


def optimize(cjson: dict, options: dict, charge: int, spin: int, debug: bool = False) -> dict:
    logger.debug("Starting geometry optimization")

    atoms_data = cjson["atoms"]
    coords_raw = atoms_data["coords"]
    if isinstance(coords_raw, dict):
        coords = coords_raw["3d"]
    else:
        coords = coords_raw
    elem = atoms_data["elements"]["number"]

    n_atoms = len(elem)
    charge_val = int(cjson.get("properties", {}).get("totalCharge", charge))
    spin_val = int(cjson.get("properties", {}).get("totalSpinMultiplicity", spin))

    mol_name = cjson.get("name", cjson.get("molecule", {}).get("name", "molecule"))
    safe_name = _sanitize_name(mol_name)

    counter = 1
    while (CALCS_DIR / f"{safe_name}_{counter:03d}").exists():
        counter += 1
    calc_dir = CALCS_DIR / f"{safe_name}_{counter:03d}"
    calc_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Calculation directory: {calc_dir}")

    xyz_lines = [f"{n_atoms}", f"{mol_name}"]
    for i in range(n_atoms):
        sym = _get_elem_symbol(elem[i])
        xyz_lines.append(
            f"{sym:3s}  {coords[3 * i]:12.8f}  {coords[3 * i + 1]:12.8f}  {coords[3 * i + 2]:12.8f}"
        )
    xyz_text = "\n".join(xyz_lines) + "\n"
    (calc_dir / "input.xyz").write_text(xyz_text, encoding="utf-8")
    logger.debug(f"Saved input.xyz ({n_atoms} atoms)")

    geom_lines = "\n".join(
        f"  {elem[i]:3d}  {coords[3 * i]:12.8f}  {coords[3 * i + 1]:12.8f}  {coords[3 * i + 2]:12.8f}"
        for i in range(n_atoms)
    )
    mol_spec = (
        f"{charge_val} {spin_val}\n{geom_lines}\nunits angstrom\nno_com\nno_reorient"
    )
    mol = psi4.geometry(mol_spec)
    mol.reset_point_group("c1")
    logger.debug(f"Psi4 molecule created (charge={charge_val}, mult={spin_val})")

    psi4.set_output_file(str(calc_dir / "psi4.log"))
    psi4_logger = logging.getLogger("psi4")
    for h in psi4_logger.handlers[:]:
        psi4_logger.removeHandler(h)
    psi4_logger.addHandler(logging.FileHandler(str(calc_dir / "psi4.log"), encoding="utf-8"))
    psi4_logger.setLevel(logging.WARNING)

    from .config import load_config
    _cfg = load_config()
    method = _option(options, "method", _cfg.get("method", "wB97X-D"))
    basis = _option(options, "basis", _cfg.get("basis", "def2-TZVP"))
    geom_maxiter = _option(options, "geom_maxiter", _cfg.get("geom_maxiter", 100))

    reference = "uhf" if (charge_val != 0 or spin_val != 1) else "rhf"
    logger.debug(f"Method: {method}, Basis: {basis}, Reference: {reference}")

    psi4.set_options({
        "basis": basis,
        "scf_type": "df",
        "reference": reference,
        "geom_maxiter": geom_maxiter,
        "e_convergence": 1e-8,
        "d_convergence": 1e-8,
    })

    converged = False
    final_energy = None
    try:
        final_energy = psi4.optimize(method, molecule=mol)
        converged = True
        logger.debug(f"Optimization converged. Final energy: {final_energy:.8f} Eh")
    except Exception as e:
        logger.warning(f"Optimization did not fully converge: {e}")
        try:
            final_energy = psi4.variable("CURRENT ENERGY")
        except Exception:
            final_energy = None
        if final_energy is None or final_energy == 0.0:
            try:
                final_energy = psi4.variable("OPTIMIZATION ENERGY")
            except Exception:
                pass

    try:
        geom = mol.geometry()
        n = geom.shape[0]
        optimized = []
        for i in range(n):
            for j in range(3):
                optimized.append(geom.get(i, j) * BOHR_TO_ANGSTROM)
        if len(optimized) != n_atoms * 3:
            raise ValueError(f"Expected {n_atoms * 3} coords, got {len(optimized)}")
        logger.debug("Extracted optimized coords from geometry matrix")
    except Exception as e:
        logger.warning(f"Could not extract optimized coordinates: {e}")
        optimized = list(coords)

    for key in ("vibrations", "basisSet", "orbitals", "cube"):
        cjson.pop(key, None)
    cjson["atoms"].pop("formalCharges", None)
    cjson["atoms"].pop("partialCharges", None)
    cjson["atoms"]["coords"] = {"3d": optimized}
    if final_energy is not None:
        cjson["properties"] = cjson.get("properties", {})
        cjson["properties"]["totalEnergy"] = round(final_energy, 8)
        cjson["properties"]["totalCharge"] = charge_val
        cjson["properties"]["totalSpinMultiplicity"] = spin_val

    opt_xyz_lines = [f"{n_atoms}", f"{mol_name} (optimized)"]
    for i in range(n_atoms):
        sym = _get_elem_symbol(elem[i])
        opt_xyz_lines.append(
            f"{sym:3s}  {optimized[3 * i]:12.8f}  {optimized[3 * i + 1]:12.8f}  {optimized[3 * i + 2]:12.8f}"
        )
    (calc_dir / "optimized.xyz").write_text("\n".join(opt_xyz_lines) + "\n", encoding="utf-8")

    energy_str = f"{final_energy:.8f} Eh" if final_energy is not None else "N/A"
    method_str = f"{method}/{basis}"
    if converged:
        message = (
            f"Geometry optimization converged.\n"
            f"Method: {method_str}\n"
            f"Final energy: {energy_str}\n"
            f"Calculations saved in: {calc_dir.name}/"
        )
    else:
        message = (
            f"Geometry optimization did NOT fully converge.\n"
            f"Method: {method_str}\n"
            f"Final energy: {energy_str}\n"
            f"Partial results saved in: {calc_dir.name}/"
        )

    logger.debug("Optimization complete, returning result")
    return {
        "moleculeFormat": "cjson",
        "cjson": cjson,
        "calcDir": str(calc_dir),
        "message": message,
    }

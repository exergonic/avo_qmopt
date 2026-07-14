import logging
import re
from contextlib import redirect_stdout

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

    psi4.set_memory("2 GB", quiet=True)
    psi4.set_output_file(str(calc_dir / "psi4.log"))
    psi4_logger = logging.getLogger("psi4")
    for h in psi4_logger.handlers[:]:
        psi4_logger.removeHandler(h)
    psi4_logger.addHandler(logging.FileHandler(str(calc_dir / "psi4.log"), encoding="utf-8"))
    psi4_logger.setLevel(logging.INFO)
    psi4_logger.propagate = False

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
        "g_convergence": "gau",
        "e_convergence": 1e-8,
        "d_convergence": 1e-8,
    })

    log_handle = open(calc_dir / "psi4.log", "a", encoding="utf-8")

    converged = False
    final_energy = None
    try:
        with redirect_stdout(log_handle):
            final_energy = psi4.optimize(method, molecule=mol)
        converged = True
        logger.debug(f"Optimization converged. Final energy: {final_energy:.8f} Eh")
    except psi4.OptimizationConvergenceError as e:
        logger.warning(f"Optimization failed to converge: {e}")
        try:
            final_energy = psi4.variable("CURRENT ENERGY")
        except Exception:
            final_energy = None
            logger.debug("Could not retrieve partial energy")
    except Exception as e:
        logger.warning(f"Optimization error: {e}")
        try:
            final_energy = psi4.variable("CURRENT ENERGY")
        except Exception:
            final_energy = None

    try:
        with redirect_stdout(log_handle):
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

    imag_msg = ""
    if _cfg.get("hess", False):
        logger.debug("Running frequency calculation after optimization")
        try:
            with redirect_stdout(log_handle):
                freq_energy, wfn = psi4.frequency(method, molecule=mol, return_wfn=True)
            fa = wfn.frequency_analysis

            omega = fa["omega"].data
            ir_raw = fa["IR_intensity"].data
            disp_mat = fa["x"].data

            n_tr = sum(1 for f in omega if abs(f) < 5.0)
            n_vib = len(omega) - n_tr
            idx_by_mag = sorted(range(len(omega)), key=lambda i: abs(omega[i]), reverse=True)
            keep = sorted(idx_by_mag[:n_vib])

            def _freq_val(f):
                return -abs(f.imag) if abs(f.imag) > abs(f.real) else f.real
            freq_list = [round(_freq_val(omega[i]), 2) for i in keep]
            ir_list = [round(float(ir_raw[i]), 4) for i in keep]

            modes_list = list(range(1, len(freq_list) + 1))
            eigen = []
            try:
                for i in keep:
                    eigen.append([float(disp_mat[i, k]) for k in range(n_atoms * 3)])
            except Exception:
                logger.debug("Normal mode eigenvectors not available from wavefunction")

            cjson["vibrations"] = {
                "frequencies": freq_list,
                "modes": modes_list,
                "intensities": ir_list,
                "eigenVectors": eigen,
            }

            n_imag = sum(1 for f in freq_list if f < -5.0)
            if n_imag > 0:
                label = "frequency" if n_imag == 1 else "frequencies"
                imag_msg = (
                    f"\n{n_imag} imaginary {label} found. "
                    "This structure may not be a true minimum on the potential energy surface."
                )
            else:
                imag_msg = "\nNo imaginary frequencies detected."

            logger.debug(f"Frequency calculation complete ({len(freq_list)} modes, {n_imag} imaginary)")
        except Exception as e:
            logger.warning(f"Frequency calculation failed: {e}")
            imag_msg = f"\nFrequency calculation failed: {e}"

    energy_str = f"{final_energy:.8f} Eh" if final_energy is not None else "N/A"
    method_str = f"{method}/{basis}"
    opt_xyz_lines = [f"{n_atoms}", f"{mol_name}  |  {method_str}  |  {energy_str}"]
    for i in range(n_atoms):
        sym = _get_elem_symbol(elem[i])
        opt_xyz_lines.append(
            f"{sym:3s}  {optimized[3 * i]:12.8f}  {optimized[3 * i + 1]:12.8f}  {optimized[3 * i + 2]:12.8f}"
        )
    (calc_dir / "optimized.xyz").write_text("\n".join(opt_xyz_lines) + "\n", encoding="utf-8")
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
    if imag_msg:
        message += imag_msg

    log_handle.close()
    logger.debug("Optimization complete, returning result")
    return {
        "moleculeFormat": "cjson",
        "cjson": cjson,
        "calcDir": str(calc_dir),
        "message": message,
    }

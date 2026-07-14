import json
from . import CALCS_DIR

CONFIG_PATH = CALCS_DIR / "config.json"

_DEFAULT_CONFIG = {
    "method": "wB97X-D",
    "basis": "def2-TZVP",
    "geom_maxiter": 100,
    "hess": False,
}

METHODS = ["HF", "B3LYP", "PBE", "PBE0", "wB97X-D"]

BASIS_SETS = [
    "cc-pVDZ",
    "aug-cc-pVDZ",
    "cc-pVTZ",
    "aug-cc-pVTZ",
    "def2-SVP",
    "def2-SVPD",
    "def2-TZVP",
    "def2-TZVPD",
]


def load_config():
    if not CONFIG_PATH.exists():
        return dict(_DEFAULT_CONFIG)
    try:
        return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return dict(_DEFAULT_CONFIG)


def save_config(config):
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CONFIG_PATH.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


def get_config_options():
    config = load_config()
    method_default = config.get("method", "HF")
    basis_default = config.get("basis", "def2-TZVP")
    return {
        "method": {
            "type": "stringList",
            "label": "Method",
            "values": METHODS,
            "default": METHODS.index(method_default)
            if method_default in METHODS
            else METHODS.index("wB97X-D"),
            "order": 1.0,
        },
        "basis": {
            "type": "stringList",
            "label": "Basis Set",
            "values": BASIS_SETS,
            "default": BASIS_SETS.index(basis_default)
            if basis_default in BASIS_SETS
            else BASIS_SETS.index("def2-TZVP"),
            "order": 2.0,
        },
        "geom_maxiter": {
            "type": "integer",
            "label": "Max Optimization Cycles",
            "minimum": 1,
            "maximum": 500,
            "default": config.get("geom_maxiter", 100),
            "order": 3.0,
        },
        "hess": {
            "type": "boolean",
            "label": "Compute vibrational frequencies after optimization",
            "default": config.get("hess", False),
            "order": 4.0,
        },
        "note": {
            "type": "text",
            "label": "Note",
            "default": (
                "\nRecommended presets:\n"
                "  General use (recommended):\t\twB97X-D / def2-TZVP\n"
                "  Fast / large systems:\t\t\tHF / def2-SVP\n"
                "  High accuracy:\t\t\twB97X-D / aug-cc-pVTZ\n"
                "  Charged / anions:\t\t\twB97X-D / def2-TZVPD\n"
                "\n"
                "Memory requirements increase with system size and basis set.\n"
                "Larger systems (30+ atoms) with triple-zeta or\n"
                "diffuse basis may require significant memory.\n"
                "Switch to a smaller basis if the calculation fails."
            ),
            "order": 99.0,
        },
    }


def update_config(avo_input):
    options = avo_input.get("options", {})
    config = load_config()
    changed = False
    for key in ("method", "basis", "geom_maxiter", "hess"):
        if key in options:
            config[key] = options[key]
            changed = True
    if changed:
        save_config(config)
    return {"cjson": avo_input.get("cjson", {})}

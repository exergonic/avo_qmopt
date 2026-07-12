import logging
import platform
import subprocess

from . import CALCS_DIR

logger = logging.getLogger(__name__)


def open_calcs_dir(cjson: dict) -> dict:
    calc_dir = CALCS_DIR.resolve()
    logger.debug(f"Opening calculations directory: {calc_dir}")
    if platform.system() == "Windows":
        subprocess.run(["explorer.exe", str(calc_dir)])
    elif platform.system() == "Darwin":
        subprocess.run(["open", str(calc_dir)])
    else:
        subprocess.run(["xdg-open", str(calc_dir)])

    return {
        "moleculeFormat": "cjson",
        "cjson": cjson,
    }

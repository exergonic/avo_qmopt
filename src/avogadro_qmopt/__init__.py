import argparse
import json
import logging
import sys
import traceback
from pathlib import Path

# Make sure stdout stream is always Unicode, as Avogadro expects
sys.stdout.reconfigure(encoding="utf-8")

logger = logging.getLogger(__name__)

def main():
    logging.basicConfig(
        level=logging.DEBUG,
        stream=sys.stderr,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    logger.debug("CLI input: " + " ".join(sys.argv))

    parser = argparse.ArgumentParser(
        description="Avogadro QM Geometry Optimization plugin"
    )
    parser.add_argument("feature", nargs="?", default="qmopt", help="Feature to run")
    parser.add_argument("--lang", default="en_US", help="Language")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    parser.add_argument("--user-options", action="store_true",
                        help="Return available config options (Avogadro protocol)")
    args = parser.parse_args()
    logger.debug(f"Parsed args: {args}")

    raw = sys.stdin.read()
    logger.debug(f"Read {len(raw)} bytes from stdin")
    data = json.loads(raw)
    logger.debug(f"Input keys: {', '.join(data.keys())}")

    cjson = data.get("cjson", {})
    options = data.get("options", {})
#    charge = data.get("charge", 0)
#    spin = data.get("spin", 1)

    try:
        if args.feature == "Test":
            result = { "message": "Test Success!" }
#       if args.feature == "config":
#           from .config import get_config_options, update_config
#           if args.user_options:
#               result = {"userOptions": get_config_options()}
#           else:
#               result = update_config(data)
 #       elif args.feature == "ibo":
 #           from .calcs import compute_ibo
 #           from .config import load_config
 #           _cfg = load_config()
 #           charge = _cfg.get("charge", charge)
 #           spin = _cfg.get("mult", spin)

 #           result = compute_ibo(cjson, options, charge, spin, debug=args.debug)
 #       elif args.feature == "open":
 #           from .links import open_calcs_dir

 #           result = open_calcs_dir(cjson)

        else:
            result = {"error": f"Unknown feature: {args.feature}"}
    except Exception as e:
        limit = None if args.debug else 3
        result = {"error": "".join(traceback.format_exception(e, limit=limit))}
        logger.exception("Unhandled exception")

    print(json.dumps(result))


if __name__ == "__main__":
    main()

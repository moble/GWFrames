"""Choose or set the package version in pyproject.toml.

    python bump_version.py              # print the version after a bump
    python bump_version.py --set 2.1.0  # write that version to pyproject.toml

The bump is chosen from the commit messages in the JSON list in the environment
variable COMMIT_MESSAGES.  The largest one mentioned in any commit wins:

    #major           2.3.4 -> 3.0.0
    #minor           2.3.4 -> 2.4.0
    #patch           2.3.4 -> 2.3.5   (the default, if nothing else is given)

Adding #dev (or #pre) makes a development release instead:

    #dev             2.3.4 -> 2.3.5.dev0,  2.3.5.dev0 -> 2.3.5.dev1
    #minor #dev      2.3.4 -> 2.4.0.dev0

and a plain #patch (or no tag at all) after a development release completes it:
2.3.5.dev1 -> 2.3.5.  When run in GitHub Actions, the new version and whether it
is a prerelease are also written to GITHUB_OUTPUT.

"""

import json
import os
import re
import sys
from pathlib import Path

PYPROJECT = Path(__file__).resolve().parents[2] / "pyproject.toml"
VERSION_LINE = re.compile(r'^version = "([^"]+)"$', re.MULTILINE)
VERSION = re.compile(r"^(\d+)\.(\d+)\.(\d+)(?:\.dev(\d+))?$")


def current_version():
    match = VERSION_LINE.search(PYPROJECT.read_text())
    if not match:
        sys.exit(f"Could not find the version in {PYPROJECT}")
    return match.group(1)


def bump(version, messages):
    match = VERSION.match(version)
    if not match:
        sys.exit(f"Cannot bump version '{version}'; expected X.Y.Z or X.Y.Z.devN")
    major, minor, patch = (int(x) for x in match.group(1, 2, 3))
    dev = match.group(4)
    text = "\n".join(messages)

    def tagged(*names):
        return any(re.search(rf"#{name}\b", text, re.IGNORECASE) for name in names)

    want_dev = tagged("dev", "pre", "prerelease")
    if tagged("major"):
        new = f"{major + 1}.0.0"
    elif tagged("minor"):
        new = f"{major}.{minor + 1}.0"
    elif dev is not None:
        # Continue the development series, or complete it
        if want_dev:
            return f"{major}.{minor}.{patch}.dev{int(dev) + 1}"
        return f"{major}.{minor}.{patch}"
    else:
        new = f"{major}.{minor}.{patch + 1}"
    return f"{new}.dev0" if want_dev else new


def write_outputs(version):
    output = os.environ.get("GITHUB_OUTPUT")
    if output:
        with open(output, "a") as f:
            f.write(f"version={version}\n")
            f.write(f"prerelease={'true' if '.dev' in version else 'false'}\n")


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--set":
        version = sys.argv[2]
        PYPROJECT.write_text(VERSION_LINE.sub(f'version = "{version}"', PYPROJECT.read_text(), count=1))
        print(version)
        return
    messages = json.loads(os.environ.get("COMMIT_MESSAGES") or "[]")
    version = bump(current_version(), messages)
    write_outputs(version)
    print(version)


if __name__ == "__main__":
    main()

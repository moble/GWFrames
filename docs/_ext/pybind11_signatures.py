"""Tidy the signatures that pybind11 writes into docstrings.

pybind11 starts each docstring with the function's signature, which autodoc
reads, but its signatures name the private extension module and the
``typing.Supports*`` protocols, and include ``self`` for methods.  For an
overloaded function, pybind11 writes ``(*args, **kwargs)`` followed by a
numbered list of the real signatures, which autodoc does not understand.  This
extension rewrites both forms, giving each overload its own signature line.
"""

import re

OVERLOADED = "Overloaded function."
OVERLOAD_LINE = re.compile(r"^\d+\. (?P<name>\w+)(?P<sig>\(.*\))(?: -> (?P<ret>.*))?$")


def simplify(text):
    """Replace pybind11's spelling of types with the public names."""
    text = text.replace("GWFrames._GWFrames.", "GWFrames.")
    text = re.sub(r"typing\.Supports(Int|Float|Complex) \| typing\.SupportsIndex", lambda m: m[1].lower(), text)
    text = re.sub(r"typing\.Supports(Int|Float|Complex)", lambda m: m[1].lower(), text)
    # The default for lists of ell values is an empty array, meaning all of them
    text = re.sub(r"array\(\[\], dtype=\w+\)", "[]", text)
    return text


def drop_self(signature):
    """Remove the leading ``self`` argument from a method signature."""
    return re.sub(r"^\(self(?:: [\w.]+)?(?:, |(?=\)))", "(", signature)


def overloads(obj):
    """Return the (signature, return annotation) pairs of an overloaded function."""
    doc = getattr(obj, "__doc__", None) or ""
    if OVERLOADED not in doc:
        return []
    return [(m["sig"], m["ret"]) for m in map(OVERLOAD_LINE.match, doc.splitlines()) if m]


def process_signature(app, what, name, obj, options, signature, return_annotation):
    if what not in ("function", "method") or signature is None:
        return None
    pairs = overloads(obj)
    if not pairs:
        if what == "method":
            signature = drop_self(signature)
        return simplify(signature), simplify(return_annotation) if return_annotation else return_annotation
    # autodoc accepts only one signature from this event, but it writes the
    # signature into the directive verbatim, so the other overloads can follow
    # on continuation lines, each repeating the name as autodoc writes it:
    # `Class.method` for a method, and the bare name for a function.
    short_name = ".".join(name.split(".")[-2 if what == "method" else -1 :])
    lines = []
    for sig, ret in pairs:
        if what == "method":
            sig = drop_self(sig)
        lines.append(simplify(sig + (" -> " + ret if ret else "")))
    return f"\n   {short_name}".join(lines), None


def process_docstring(app, what, name, obj, options, lines):
    """Remove pybind11's list of overloads, now shown as signatures."""
    if OVERLOADED not in lines:
        return
    kept = [line for line in lines if line != OVERLOADED and not OVERLOAD_LINE.match(line)]
    while kept and not kept[0].strip():
        kept.pop(0)
    lines[:] = kept


def resolve_relative_to_package(app, env, node, contnode):
    """Resolve references like `Quaternions.Squad` relative to the GWFrames package.

    The docstrings name objects as they are written in Python after
    ``import GWFrames``, which does not resolve from within the
    `GWFrames.Quaternions` submodule.
    """
    if node.get("refdomain") != "py" or node["reftarget"].startswith("GWFrames."):
        return None
    domain = env.domains["py"]
    target = node["reftarget"]
    node["reftarget"] = "GWFrames." + target
    try:
        return domain.resolve_xref(env, node["refdoc"], app.builder, node["reftype"], node["reftarget"], node, contnode)
    finally:
        node["reftarget"] = target


def setup(app):
    app.connect("autodoc-process-signature", process_signature)
    app.connect("autodoc-process-docstring", process_docstring)
    app.connect("missing-reference", resolve_relative_to_package)
    return {"parallel_read_safe": True, "parallel_write_safe": True}

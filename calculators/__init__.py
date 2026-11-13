"""
automag.calculators
===================

Registry of calculator interfaces.

The pipeline never checks ``if calculator == "vasp"``.  It asks the registry
for an interface by name and calls its methods.  Adding support for a new code
therefore requires **no changes to the pipeline or stage code**:

* either point the input file at any ASE calculator class with
  ``calculator = "ase:package.module.Calculator"`` (see
  :class:`automag.calculators.generic.GenericASEInterface`);
* or write a plugin module defining a subclass of
  :class:`~automag.calculators.base.CalculatorInterface`, decorate it with
  :func:`register`, and list the module in ``calculator_plugins`` in the input
  file (or in the ``AUTOMAG_CALCULATORS`` environment variable, comma
  separated).

.. codeauthor:: Automag contributors
"""

import importlib
import os

from calculators.base import CalculatorInterface

_REGISTRY = {}
_BUILTINS_LOADED = False


def register(cls=None):
    """Register a calculator interface class (usable as a decorator)."""

    def _register(interface):
        if not getattr(interface, "name", None):
            raise ValueError("A calculator interface must define a 'name'.")
        _REGISTRY[interface.name] = interface
        return interface

    if cls is None:
        return _register
    return _register(cls)


def load_plugins(extra=None):
    """Load built-in interfaces and any external plugin modules."""
    global _BUILTINS_LOADED
    if not _BUILTINS_LOADED:
        from calculators.vasp import VaspInterface
        from calculators.espresso import EspressoInterface
        from calculators.aims import FhiAimsInterface
        from calculators.prophet import ProphetSpinInterface
        from calculators.generic import GenericASEInterface

        _REGISTRY.setdefault("vasp", VaspInterface)
        _REGISTRY.setdefault("espresso", EspressoInterface)
        _REGISTRY.setdefault("qe", EspressoInterface)
        _REGISTRY.setdefault("fhi-aims", FhiAimsInterface)
        _REGISTRY.setdefault("aims", FhiAimsInterface)
        _REGISTRY.setdefault("prophet_spin", ProphetSpinInterface)
        _REGISTRY.setdefault("prophet", ProphetSpinInterface)
        _REGISTRY.setdefault("ase", GenericASEInterface)
        _BUILTINS_LOADED = True

    modules = list(extra or [])
    env = os.environ.get("AUTOMAG_CALCULATORS", "")
    modules.extend(item.strip() for item in env.split(",") if item.strip())
    for module_name in modules:
        importlib.import_module(module_name)


def get_calculator(name, **kwargs):
    """Return an instance of the interface registered under ``name``."""
    load_plugins(kwargs.pop("plugins", None))
    if name in _REGISTRY:
        return _REGISTRY[name](**kwargs)
    if name.startswith("ase:"):
        from calculators.generic import GenericASEInterface

        return GenericASEInterface(calculator_path=name[4:], **kwargs)
    available = ", ".join(sorted(_REGISTRY))
    raise KeyError(
        f"Unknown calculator '{name}'. Available: {available}. "
        "Use 'ase:<dotted.path>' for any ASE calculator or register a plugin."
    )


def available():
    """Return the names of all registered interfaces."""
    load_plugins()
    return sorted(_REGISTRY)

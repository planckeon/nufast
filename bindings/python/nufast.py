"""
NuFast Python Bindings
======================

Fast and accurate three-flavor neutrino oscillation probabilities in vacuum and matter.

This module provides Python bindings to the NuFast Zig library via ctypes.
It implements the algorithm from Denton & Parke (arXiv:2405.02400).

Example
-------
>>> import nufast
>>> 
>>> # Vacuum oscillation at DUNE baseline
>>> probs = nufast.vacuum_probability(L=1300, E=2.5)
>>> print(f"P(νμ → νe) = {probs.Pme:.4f}")
>>> 
>>> # Matter effects
>>> probs = nufast.matter_probability(L=1300, E=2.5, rho=2.848)
>>> print(f"P(νμ → νe) with matter = {probs.Pme:.4f}")

"""

from __future__ import annotations

import ctypes
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

# =============================================================================
# Library Loading
# =============================================================================

def _find_library() -> Path:
    """Find the nufast shared library."""
    # Check for explicit path via environment variable
    if lib_path := os.environ.get("NUFAST_LIB"):
        return Path(lib_path)
    
    # Platform-specific library names
    if sys.platform == "win32":
        lib_names = ["nufast.dll", "libnufast.dll"]
    elif sys.platform == "darwin":
        lib_names = ["libnufast.dylib", "libnufast.so"]
    else:
        lib_names = ["libnufast.so"]
    
    # Search paths (relative to this file)
    this_dir = Path(__file__).parent
    search_dirs = [
        this_dir,
        this_dir / "lib",
        this_dir.parent.parent / "benchmarks" / "zig" / "zig-out" / "lib",
        Path.cwd(),
        Path.cwd() / "lib",
    ]
    
    for search_dir in search_dirs:
        for lib_name in lib_names:
            lib_path = search_dir / lib_name
            if lib_path.exists():
                return lib_path
    
    raise FileNotFoundError(
        f"Could not find nufast library. Searched for {lib_names} in {search_dirs}. "
        "Set NUFAST_LIB environment variable to the library path, or run 'zig build lib' "
        "in the benchmarks/zig directory."
    )


def _load_library() -> ctypes.CDLL:
    """Load the nufast shared library."""
    lib_path = _find_library()
    lib = ctypes.CDLL(str(lib_path))
    
    # =========================================================================
    # Function signatures
    # =========================================================================
    
    # Parameter setters
    lib.nufast_set_vacuum_params.argtypes = [
        ctypes.c_double,  # s12sq
        ctypes.c_double,  # s13sq
        ctypes.c_double,  # s23sq
        ctypes.c_double,  # delta
        ctypes.c_double,  # Dmsq21
        ctypes.c_double,  # Dmsq31
        ctypes.c_bool,    # antineutrino
    ]
    lib.nufast_set_vacuum_params.restype = None
    
    lib.nufast_set_default_params.argtypes = []
    lib.nufast_set_default_params.restype = None
    
    lib.nufast_set_matter_params.argtypes = [
        ctypes.c_double,  # rho
        ctypes.c_double,  # Ye
        ctypes.c_uint8,   # n_newton
        ctypes.c_bool,    # antineutrino
    ]
    lib.nufast_set_matter_params.restype = None
    
    # Probability calculations (stateful)
    lib.nufast_vacuum_probability.argtypes = [ctypes.c_double, ctypes.c_double]
    lib.nufast_vacuum_probability.restype = ctypes.POINTER(ctypes.c_double)
    
    lib.nufast_matter_probability.argtypes = [ctypes.c_double, ctypes.c_double]
    lib.nufast_matter_probability.restype = ctypes.POINTER(ctypes.c_double)
    
    # Direct calculations (stateless)
    lib.nufast_vacuum_prob_direct.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,  # mixing
        ctypes.c_double, ctypes.c_double,  # mass differences
        ctypes.c_double, ctypes.c_double,  # L, E
        ctypes.c_bool,    # antineutrino
        ctypes.POINTER(ctypes.c_double),  # output buffer
    ]
    lib.nufast_vacuum_prob_direct.restype = None
    
    lib.nufast_matter_prob_direct.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,  # mixing
        ctypes.c_double, ctypes.c_double,  # mass differences
        ctypes.c_double, ctypes.c_double,  # L, E
        ctypes.c_double, ctypes.c_double,  # rho, Ye
        ctypes.c_uint8,   # n_newton
        ctypes.c_bool,    # antineutrino
        ctypes.POINTER(ctypes.c_double),  # output buffer
    ]
    lib.nufast_matter_prob_direct.restype = None
    
    # Quick Pme calculations
    lib.nufast_vacuum_Pme.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double,
    ]
    lib.nufast_vacuum_Pme.restype = ctypes.c_double
    
    lib.nufast_vacuum_Pme_default.argtypes = [ctypes.c_double, ctypes.c_double]
    lib.nufast_vacuum_Pme_default.restype = ctypes.c_double
    
    lib.nufast_matter_Pme_default.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
    lib.nufast_matter_Pme_default.restype = ctypes.c_double
    
    # Batch processing
    lib.nufast_init_vacuum_batch.argtypes = []
    lib.nufast_init_vacuum_batch.restype = None
    
    lib.nufast_init_matter_batch.argtypes = []
    lib.nufast_init_matter_batch.restype = None
    
    lib.nufast_get_energies_ptr.argtypes = []
    lib.nufast_get_energies_ptr.restype = ctypes.POINTER(ctypes.c_double)
    
    lib.nufast_get_batch_output_ptr.argtypes = []
    lib.nufast_get_batch_output_ptr.restype = ctypes.POINTER(ctypes.c_double)
    
    lib.nufast_get_max_batch_size.argtypes = []
    lib.nufast_get_max_batch_size.restype = ctypes.c_size_t
    
    lib.nufast_vacuum_batch_Pme.argtypes = [ctypes.c_double, ctypes.c_size_t]
    lib.nufast_vacuum_batch_Pme.restype = ctypes.POINTER(ctypes.c_double)
    
    lib.nufast_matter_batch_Pme.argtypes = [ctypes.c_double, ctypes.c_size_t]
    lib.nufast_matter_batch_Pme.restype = ctypes.POINTER(ctypes.c_double)
    
    # Stateless batch
    lib.nufast_vacuum_batch_Pme_direct.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double,
        ctypes.c_double,  # L
        ctypes.POINTER(ctypes.c_double),  # energies
        ctypes.c_size_t,  # count
        ctypes.c_bool,    # antineutrino
        ctypes.POINTER(ctypes.c_double),  # output
    ]
    lib.nufast_vacuum_batch_Pme_direct.restype = None
    
    lib.nufast_matter_batch_Pme_direct.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double,
        ctypes.c_double,  # L
        ctypes.POINTER(ctypes.c_double),  # energies
        ctypes.c_size_t,  # count
        ctypes.c_double, ctypes.c_double,  # rho, Ye
        ctypes.c_uint8,   # n_newton
        ctypes.c_bool,    # antineutrino
        ctypes.POINTER(ctypes.c_double),  # output
    ]
    lib.nufast_matter_batch_Pme_direct.restype = None
    
    # Default parameter accessors
    lib.nufast_default_s12sq.argtypes = []
    lib.nufast_default_s12sq.restype = ctypes.c_double
    
    lib.nufast_default_s13sq.argtypes = []
    lib.nufast_default_s13sq.restype = ctypes.c_double
    
    lib.nufast_default_s23sq.argtypes = []
    lib.nufast_default_s23sq.restype = ctypes.c_double
    
    lib.nufast_default_delta.argtypes = []
    lib.nufast_default_delta.restype = ctypes.c_double
    
    lib.nufast_default_Dmsq21.argtypes = []
    lib.nufast_default_Dmsq21.restype = ctypes.c_double
    
    lib.nufast_default_Dmsq31.argtypes = []
    lib.nufast_default_Dmsq31.restype = ctypes.c_double
    
    lib.nufast_default_rho.argtypes = []
    lib.nufast_default_rho.restype = ctypes.c_double
    
    lib.nufast_default_Ye.argtypes = []
    lib.nufast_default_Ye.restype = ctypes.c_double
    
    return lib


# Load the library lazily
_lib: ctypes.CDLL | None = None


def _get_lib() -> ctypes.CDLL:
    """Get the loaded library, loading it if necessary."""
    global _lib
    if _lib is None:
        _lib = _load_library()
    return _lib


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class OscillationParams:
    """Neutrino oscillation parameters.
    
    Attributes
    ----------
    s12sq : float
        sin²θ₁₂ (solar mixing angle squared)
    s13sq : float
        sin²θ₁₃ (reactor mixing angle squared)
    s23sq : float
        sin²θ₂₃ (atmospheric mixing angle squared)
    delta : float
        CP-violating phase δ (radians)
    Dmsq21 : float
        Δm²₂₁ (solar mass-squared difference, eV²)
    Dmsq31 : float
        Δm²₃₁ (atmospheric mass-squared difference, eV²)
        Positive for normal ordering, negative for inverted ordering.
    """
    s12sq: float
    s13sq: float
    s23sq: float
    delta: float
    Dmsq21: float
    Dmsq31: float
    
    @classmethod
    def default(cls) -> OscillationParams:
        """Return default NuFIT 5.2 parameters (normal ordering)."""
        lib = _get_lib()
        return cls(
            s12sq=lib.nufast_default_s12sq(),
            s13sq=lib.nufast_default_s13sq(),
            s23sq=lib.nufast_default_s23sq(),
            delta=lib.nufast_default_delta(),
            Dmsq21=lib.nufast_default_Dmsq21(),
            Dmsq31=lib.nufast_default_Dmsq31(),
        )


@dataclass
class ProbabilityMatrix:
    """3×3 neutrino oscillation probability matrix.
    
    Elements are P(α → β) where α is the initial flavor and β is the final flavor.
    Rows are indexed by initial flavor (e, μ, τ), columns by final flavor.
    
    Attributes
    ----------
    Pee, Pem, Pet : float
        Electron neutrino row
    Pme, Pmm, Pmt : float
        Muon neutrino row
    Pte, Ptm, Ptt : float
        Tau neutrino row
    """
    Pee: float
    Pem: float
    Pet: float
    Pme: float
    Pmm: float
    Pmt: float
    Pte: float
    Ptm: float
    Ptt: float
    
    def __getitem__(self, key: tuple[int, int]) -> float:
        """Get element by index: matrix[i, j] = P(flavor_i → flavor_j)."""
        i, j = key
        return self.as_tuple()[i][j]
    
    def as_tuple(self) -> tuple[tuple[float, float, float], ...]:
        """Return as nested tuple ((Pee, Pem, Pet), (Pme, Pmm, Pmt), (Pte, Ptm, Ptt))."""
        return (
            (self.Pee, self.Pem, self.Pet),
            (self.Pme, self.Pmm, self.Pmt),
            (self.Pte, self.Ptm, self.Ptt),
        )
    
    def as_list(self) -> list[list[float]]:
        """Return as nested list for numpy compatibility."""
        return [list(row) for row in self.as_tuple()]
    
    @classmethod
    def from_buffer(cls, buf: ctypes.POINTER) -> ProbabilityMatrix:
        """Create from a ctypes buffer (9 f64s, row-major)."""
        return cls(
            Pee=buf[0], Pem=buf[1], Pet=buf[2],
            Pme=buf[3], Pmm=buf[4], Pmt=buf[5],
            Pte=buf[6], Ptm=buf[7], Ptt=buf[8],
        )


# =============================================================================
# Experiment Presets
# =============================================================================

@dataclass
class Experiment:
    """Configuration for a neutrino oscillation experiment.
    
    Attributes
    ----------
    name : str
        Experiment name
    L : float
        Baseline distance (km)
    E : float
        Typical/peak neutrino energy (GeV)
    rho : float
        Average matter density along baseline (g/cm³)
    Ye : float
        Electron fraction (typically 0.5 for Earth)
    """
    name: str
    L: float
    E: float
    rho: float
    Ye: float = 0.5


# Common experiment configurations
EXPERIMENTS = {
    "T2K": Experiment("T2K", L=295.0, E=0.6, rho=2.6),
    "NOvA": Experiment("NOvA", L=810.0, E=2.0, rho=2.84),
    "DUNE": Experiment("DUNE", L=1300.0, E=2.5, rho=2.848),
    "Hyper-K": Experiment("Hyper-K", L=295.0, E=0.6, rho=2.6),
    "JUNO": Experiment("JUNO", L=52.5, E=0.004, rho=2.6),  # 4 MeV
}


# =============================================================================
# Main API
# =============================================================================

def vacuum_probability(
    L: float,
    E: float,
    *,
    params: OscillationParams | None = None,
    antineutrino: bool = False,
) -> ProbabilityMatrix:
    """Calculate vacuum neutrino oscillation probabilities.
    
    Parameters
    ----------
    L : float
        Baseline distance (km)
    E : float
        Neutrino energy (GeV)
    params : OscillationParams, optional
        Oscillation parameters. If None, uses NuFIT 5.2 defaults.
    antineutrino : bool, default False
        If True, calculate for antineutrinos (flips sign of δ)
    
    Returns
    -------
    ProbabilityMatrix
        3×3 oscillation probability matrix
    
    Example
    -------
    >>> probs = vacuum_probability(1300, 2.5)
    >>> print(f"P(νμ → νe) = {probs.Pme:.4f}")
    """
    lib = _get_lib()
    
    if params is None:
        params = OscillationParams.default()
    
    # Use stateless direct API
    out = (ctypes.c_double * 9)()
    lib.nufast_vacuum_prob_direct(
        params.s12sq, params.s13sq, params.s23sq, params.delta,
        params.Dmsq21, params.Dmsq31,
        L, E, antineutrino, out
    )
    
    return ProbabilityMatrix.from_buffer(out)


def matter_probability(
    L: float,
    E: float,
    *,
    rho: float | None = None,
    Ye: float = 0.5,
    n_newton: int = 0,
    params: OscillationParams | None = None,
    antineutrino: bool = False,
) -> ProbabilityMatrix:
    """Calculate matter-affected neutrino oscillation probabilities.
    
    Parameters
    ----------
    L : float
        Baseline distance (km)
    E : float
        Neutrino energy (GeV)
    rho : float, optional
        Matter density (g/cm³). Default: 2.848 (DUNE-like)
    Ye : float, default 0.5
        Electron fraction
    n_newton : int, default 0
        Newton-Raphson iterations (0-3). Higher = more accurate but slower.
    params : OscillationParams, optional
        Oscillation parameters. If None, uses NuFIT 5.2 defaults.
    antineutrino : bool, default False
        If True, calculate for antineutrinos (flips matter potential and δ)
    
    Returns
    -------
    ProbabilityMatrix
        3×3 oscillation probability matrix with matter effects
    
    Example
    -------
    >>> probs = matter_probability(1300, 2.5, rho=2.848)
    >>> print(f"P(νμ → νe) = {probs.Pme:.4f}")
    """
    lib = _get_lib()
    
    if params is None:
        params = OscillationParams.default()
    
    if rho is None:
        rho = lib.nufast_default_rho()
    
    out = (ctypes.c_double * 9)()
    lib.nufast_matter_prob_direct(
        params.s12sq, params.s13sq, params.s23sq, params.delta,
        params.Dmsq21, params.Dmsq31,
        L, E, rho, Ye, n_newton, antineutrino, out
    )
    
    return ProbabilityMatrix.from_buffer(out)


def vacuum_Pme(
    L: float,
    E: float,
    *,
    params: OscillationParams | None = None,
) -> float:
    """Quick calculation of P(νμ → νe) in vacuum.
    
    This is the most common appearance probability. For the full matrix,
    use `vacuum_probability()`.
    
    Parameters
    ----------
    L : float
        Baseline distance (km)
    E : float
        Neutrino energy (GeV)
    params : OscillationParams, optional
        Oscillation parameters. If None, uses defaults.
    
    Returns
    -------
    float
        P(νμ → νe)
    """
    lib = _get_lib()
    
    if params is None:
        return lib.nufast_vacuum_Pme_default(L, E)
    
    return lib.nufast_vacuum_Pme(
        params.s12sq, params.s13sq, params.s23sq, params.delta,
        params.Dmsq21, params.Dmsq31, L, E
    )


def matter_Pme(
    L: float,
    E: float,
    *,
    rho: float | None = None,
    Ye: float = 0.5,
    n_newton: int = 0,
    params: OscillationParams | None = None,
    antineutrino: bool = False,
) -> float:
    """Quick calculation of P(νμ → νe) with matter effects.
    
    Parameters
    ----------
    L : float
        Baseline distance (km)
    E : float
        Neutrino energy (GeV)
    rho : float, optional
        Matter density (g/cm³). Default: 2.848
    Ye : float, default 0.5
        Electron fraction
    n_newton : int, default 0
        Newton-Raphson iterations (0-3)
    params : OscillationParams, optional
        Oscillation parameters. If None, uses defaults.
    antineutrino : bool, default False
        If True, calculate for antineutrinos
    
    Returns
    -------
    float
        P(νμ → νe)
    """
    lib = _get_lib()
    
    if params is None and rho is None:
        return lib.nufast_matter_Pme_default(L, E, lib.nufast_default_rho())
    
    probs = matter_probability(L, E, rho=rho, Ye=Ye, n_newton=n_newton,
                               params=params, antineutrino=antineutrino)
    return probs.Pme


# =============================================================================
# Batch Processing
# =============================================================================

def vacuum_Pme_batch(
    L: float,
    energies: Sequence[float],
    *,
    params: OscillationParams | None = None,
    antineutrino: bool = False,
) -> list[float]:
    """Calculate P(νμ → νe) for multiple energies (vacuum).
    
    Optimized for calculating the same baseline with many energies.
    
    Parameters
    ----------
    L : float
        Baseline distance (km)
    energies : sequence of float
        Neutrino energies (GeV)
    params : OscillationParams, optional
        Oscillation parameters. If None, uses defaults.
    antineutrino : bool, default False
        If True, calculate for antineutrinos
    
    Returns
    -------
    list of float
        P(νμ → νe) for each energy
    """
    lib = _get_lib()
    
    if params is None:
        params = OscillationParams.default()
    
    count = len(energies)
    energies_arr = (ctypes.c_double * count)(*energies)
    out = (ctypes.c_double * count)()
    
    lib.nufast_vacuum_batch_Pme_direct(
        params.s12sq, params.s13sq, params.s23sq, params.delta,
        params.Dmsq21, params.Dmsq31,
        L, energies_arr, count, antineutrino, out
    )
    
    return list(out)


def matter_Pme_batch(
    L: float,
    energies: Sequence[float],
    *,
    rho: float | None = None,
    Ye: float = 0.5,
    n_newton: int = 0,
    params: OscillationParams | None = None,
    antineutrino: bool = False,
) -> list[float]:
    """Calculate P(νμ → νe) for multiple energies (with matter effects).
    
    Optimized for calculating the same baseline with many energies.
    
    Parameters
    ----------
    L : float
        Baseline distance (km)
    energies : sequence of float
        Neutrino energies (GeV)
    rho : float, optional
        Matter density (g/cm³). Default: 2.848
    Ye : float, default 0.5
        Electron fraction
    n_newton : int, default 0
        Newton-Raphson iterations (0-3)
    params : OscillationParams, optional
        Oscillation parameters. If None, uses defaults.
    antineutrino : bool, default False
        If True, calculate for antineutrinos
    
    Returns
    -------
    list of float
        P(νμ → νe) for each energy
    """
    lib = _get_lib()
    
    if params is None:
        params = OscillationParams.default()
    
    if rho is None:
        rho = lib.nufast_default_rho()
    
    count = len(energies)
    energies_arr = (ctypes.c_double * count)(*energies)
    out = (ctypes.c_double * count)()
    
    lib.nufast_matter_batch_Pme_direct(
        params.s12sq, params.s13sq, params.s23sq, params.delta,
        params.Dmsq21, params.Dmsq31,
        L, energies_arr, count, rho, Ye, n_newton, antineutrino, out
    )
    
    return list(out)


# =============================================================================
# Convenience Functions
# =============================================================================

def experiment_probability(
    experiment: Experiment | str,
    *,
    E: float | None = None,
    params: OscillationParams | None = None,
    antineutrino: bool = False,
    n_newton: int = 0,
) -> ProbabilityMatrix:
    """Calculate oscillation probabilities for a preset experiment.
    
    Parameters
    ----------
    experiment : Experiment or str
        Experiment configuration or name ("DUNE", "T2K", "NOvA", etc.)
    E : float, optional
        Override experiment's default energy (GeV)
    params : OscillationParams, optional
        Oscillation parameters. If None, uses defaults.
    antineutrino : bool, default False
        If True, calculate for antineutrinos
    n_newton : int, default 0
        Newton-Raphson iterations for matter calculation
    
    Returns
    -------
    ProbabilityMatrix
        3×3 oscillation probability matrix with matter effects
    
    Example
    -------
    >>> probs = experiment_probability("DUNE")
    >>> print(f"P(νμ → νe) at DUNE = {probs.Pme:.4f}")
    """
    if isinstance(experiment, str):
        experiment = EXPERIMENTS[experiment]
    
    energy = E if E is not None else experiment.E
    
    return matter_probability(
        experiment.L, energy,
        rho=experiment.rho, Ye=experiment.Ye,
        n_newton=n_newton, params=params, antineutrino=antineutrino
    )


# =============================================================================
# Module-level convenience
# =============================================================================

# Pre-load library defaults at import time (optional)
try:
    DEFAULT_PARAMS = OscillationParams.default()
except FileNotFoundError:
    DEFAULT_PARAMS = None  # Library not found; will error on first use

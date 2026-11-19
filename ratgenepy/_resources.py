"""Helpers for limiting CPU threads and total memory usage.

The solver backends have their own ``threads`` and ``memory`` parameters.  In
addition, :func:`parallel_limits` provides a process-wide guard used by
:func:`ratgenepy.ratio.ratioGenePar`:

* common BLAS/OpenMP/NumExpr environment variables are pinned to the requested
  thread count, so NumPy/SciPy operations do not spawn all available cores;
* when ``threadpoolctl`` is installed, already-loaded native thread pools are
  limited as well;
* on POSIX systems ``resource.setrlimit(RLIMIT_AS, ...)`` enforces an
  address-space limit for the duration of the call; on Windows a Job Object
  process-memory limit is used instead.

Memory values are expressed in megabytes (MB).  ``float("inf")`` means "do not
limit memory".
"""

from __future__ import annotations

import contextlib
import math
import os
from typing import Iterator, Tuple

DEFAULT_THREADS = 1
DEFAULT_MEMORY = float("inf")

# Variables honoured by the most common numerical libraries.
_THREAD_ENVIRONMENT_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "GOTO_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMEXPR_MAX_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
    "CRAYBLAS_NUM_THREADS",
    "TBB_NUM_THREADS",
    "RAYON_NUM_THREADS",
    "POLARS_MAX_THREADS",
)

# Prevent OpenMP/MKL from silently changing the requested thread count.
_THREAD_CONTROL_ENVIRONMENT_VARIABLES = {
    "OMP_DYNAMIC": "FALSE",
    "MKL_DYNAMIC": "FALSE",
}

_MEMORY_SUFFIXES = (
    ("gib", 1024.0),
    ("gb", 1024.0),
    ("g", 1024.0),
    ("mib", 1.0),
    ("mb", 1.0),
    ("m", 1.0),
    ("kib", 1.0 / 1024.0),
    ("kb", 1.0 / 1024.0),
    ("k", 1.0 / 1024.0),
)


def normalize_threads(threads=DEFAULT_THREADS) -> int:
    """Return a validated positive thread count."""
    if threads is None:
        threads = DEFAULT_THREADS
    try:
        value = int(threads)
    except (TypeError, ValueError) as exc:
        raise ValueError("threads must be a positive integer") from exc
    if value < 1:
        raise ValueError("threads must be a positive integer")
    return value


def normalize_memory(memory=DEFAULT_MEMORY) -> float:
    """Return a validated memory limit in MB.

    ``None`` and the strings ``"inf"``/``"infinity"``/``"unlimited"`` mean no
    limit.  Strings with ``K``/``M``/``G`` (or ``KiB``/``MiB``/``GiB``)
    suffixes are accepted and converted to MB.
    """
    if memory is None:
        return DEFAULT_MEMORY

    if isinstance(memory, str):
        text = memory.strip().lower().replace("_", "").replace(" ", "")
        if text in ("", "inf", "+inf", "infinity", "+infinity", "none",
                    "unlimited", "nolimit"):
            return DEFAULT_MEMORY
        multiplier = 1.0
        for suffix, factor in _MEMORY_SUFFIXES:
            if text.endswith(suffix):
                text = text[:-len(suffix)]
                multiplier = factor
                break
        try:
            memory = float(text) * multiplier
        except ValueError as exc:
            raise ValueError(
                "memory must be a non-negative number of megabytes, a string "
                "with a K/M/G suffix, or infinity"
            ) from exc

    try:
        value = float(memory)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "memory must be a non-negative number of megabytes or infinity"
        ) from exc

    if math.isnan(value) or value < 0:
        raise ValueError(
            "memory must be a non-negative number of megabytes or infinity"
        )
    return value


def _restore_environment(saved: dict) -> None:
    for name, value in saved.items():
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value


@contextlib.contextmanager
def thread_limits(threads=DEFAULT_THREADS) -> Iterator[int]:
    """Temporarily pin common numerical libraries to ``threads`` threads."""
    threads = normalize_threads(threads)

    saved = {
        name: os.environ.get(name)
        for name in _THREAD_ENVIRONMENT_VARIABLES
    }
    saved.update({
        name: os.environ.get(name)
        for name in _THREAD_CONTROL_ENVIRONMENT_VARIABLES
    })

    for name in _THREAD_ENVIRONMENT_VARIABLES:
        os.environ[name] = str(threads)
    os.environ.update(_THREAD_CONTROL_ENVIRONMENT_VARIABLES)

    try:
        try:
            from threadpoolctl import threadpool_limits
        except ImportError:
            threadpool_limits = None

        stack = contextlib.ExitStack()
        try:
            if threadpool_limits is not None:
                try:
                    stack.enter_context(threadpool_limits(limits=threads))
                except Exception:
                    # threadpoolctl can fail for unusual native libraries; the
                    # environment variables still provide a best-effort limit.
                    pass
            yield threads
        finally:
            stack.close()
    finally:
        _restore_environment(saved)


@contextlib.contextmanager
def _posix_memory_limit(limit_bytes: int) -> Iterator[None]:
    try:
        import resource
    except ImportError:  # pragma: no cover - only on unusual platforms
        yield
        return

    candidates = [
        getattr(resource, "RLIMIT_AS", None),
        getattr(resource, "RLIMIT_DATA", None),
    ]
    which = next((item for item in candidates if item is not None), None)
    if which is None:  # pragma: no cover - unusual platform
        yield
        return

    try:
        old_soft, old_hard = resource.getrlimit(which)
        if old_hard == resource.RLIM_INFINITY:
            new_soft = limit_bytes
        else:
            new_soft = min(limit_bytes, old_hard)
        resource.setrlimit(which, (new_soft, old_hard))
    except (OSError, ValueError):
        # The platform may not allow address-space limits; the solver-level
        # memory parameters still apply.
        yield
        return

    try:
        yield
    finally:
        try:
            resource.setrlimit(which, (old_soft, old_hard))
        except (OSError, ValueError):
            pass


@contextlib.contextmanager
def _windows_memory_limit(limit_bytes: int) -> Iterator[None]:
    """Apply a Job Object process-memory limit on Windows."""
    try:
        import ctypes
        from ctypes import wintypes
    except Exception:  # pragma: no cover - only on Windows
        yield
        return

    class _IoCounters(ctypes.Structure):
        _fields_ = [
            ("ReadOperationCount", ctypes.c_ulonglong),
            ("WriteOperationCount", ctypes.c_ulonglong),
            ("OtherOperationCount", ctypes.c_ulonglong),
            ("ReadTransferCount", ctypes.c_ulonglong),
            ("WriteTransferCount", ctypes.c_ulonglong),
            ("OtherTransferCount", ctypes.c_ulonglong),
        ]

    class _BasicLimitInformation(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_longlong),
            ("PerJobUserTimeLimit", ctypes.c_longlong),
            ("LimitFlags", ctypes.c_uint32),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", ctypes.c_uint32),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", ctypes.c_uint32),
            ("SchedulingClass", ctypes.c_uint32),
        ]

    class _ExtendedLimitInformation(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", _BasicLimitInformation),
            ("IoInfo", _IoCounters),
            ("ProcessMemoryLimit", ctypes.c_size_t),
            ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t),
            ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
    kernel32.CreateJobObjectW.restype = wintypes.HANDLE
    kernel32.SetInformationJobObject.argtypes = [
        wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, ctypes.c_uint32
    ]
    kernel32.SetInformationJobObject.restype = wintypes.BOOL
    kernel32.AssignProcessToJobObject.argtypes = [
        wintypes.HANDLE, wintypes.HANDLE
    ]
    kernel32.AssignProcessToJobObject.restype = wintypes.BOOL
    kernel32.GetCurrentProcess.restype = wintypes.HANDLE
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel32.CloseHandle.restype = wintypes.BOOL

    job = kernel32.CreateJobObjectW(None, None)
    if not job:  # pragma: no cover - Windows API failure
        yield
        return

    try:
        info = _ExtendedLimitInformation()
        # JOB_OBJECT_LIMIT_PROCESS_MEMORY | JOB_OBJECT_LIMIT_JOB_MEMORY
        info.BasicLimitInformation.LimitFlags = 0x00000100 | 0x00000200
        info.ProcessMemoryLimit = limit_bytes
        info.JobMemoryLimit = limit_bytes
        if not kernel32.SetInformationJobObject(
            job, 9, ctypes.byref(info), ctypes.sizeof(info)
        ):  # JobObjectExtendedLimitInformation == 9
            yield
            return
        if not kernel32.AssignProcessToJobObject(
            job, kernel32.GetCurrentProcess()
        ):
            # The current process may already belong to a non-breakaway job.
            # In that case the solver-level limit remains in effect.
            yield
            return
        yield
    finally:
        # Clear the limit first so a process that cannot leave the job (for
        # example because it is nested in another job) is not left capped.
        try:
            info.BasicLimitInformation.LimitFlags = 0
            info.ProcessMemoryLimit = 0
            info.JobMemoryLimit = 0
            kernel32.SetInformationJobObject(
                job, 9, ctypes.byref(info), ctypes.sizeof(info)
            )
        except Exception:
            pass
        # Closing the last handle to the job releases the process limit.
        kernel32.CloseHandle(job)


@contextlib.contextmanager
def memory_limit(memory=DEFAULT_MEMORY) -> Iterator[float]:
    """Temporarily enforce a process address-space limit (memory in MB)."""
    limit = normalize_memory(memory)
    if not math.isfinite(limit):
        yield limit
        return

    limit_bytes = max(0, int(limit * 1024.0 * 1024.0))
    if os.name == "nt":
        context = _windows_memory_limit(limit_bytes)
    else:
        context = _posix_memory_limit(limit_bytes)

    with context:
        yield limit


@contextlib.contextmanager
def parallel_limits(threads=DEFAULT_THREADS,
                    memory=DEFAULT_MEMORY) -> Iterator[Tuple[int, float]]:
    """Pin numerical-library threads and enforce a total-memory limit."""
    threads = normalize_threads(threads)
    memory = normalize_memory(memory)

    with thread_limits(threads):
        with memory_limit(memory):
            yield threads, memory


__all__ = [
    "DEFAULT_MEMORY",
    "DEFAULT_THREADS",
    "memory_limit",
    "normalize_memory",
    "normalize_threads",
    "parallel_limits",
    "thread_limits",
]

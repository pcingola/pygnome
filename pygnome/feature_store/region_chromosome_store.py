import importlib
from typing import Callable, Iterable

import numpy as np

from pygnome.feature_store.chromosome_feature_store import ChromosomeFeatureStore
from pygnome.genomics import GenomicFeature

MAX_LOOK_FORWARD = 1000

class FeatureField:
    def __init__(self, name: str, dtype=np.object_, default=None):
        self.name = name
        self.dtype = dtype
        self.default = default


class RegionFeatureCounter:
    """
    First-pass counter for feature counts and max length per bin.
    """
    def __init__(self, bin_size: int = 100_000):
        self.bin_size = bin_size
        self.feature_counts: dict[str, int] = {}
        self.max_lengths_by_bin: dict[str, dict[int, int]] = {}

    def add(self, feature) -> None:
        chrom = feature.chrom
        self.feature_counts.setdefault(chrom, 0)
        self.max_lengths_by_bin.setdefault(chrom, {})
        self.feature_counts[chrom] += 1
        bin_id = feature.start // self.bin_size
        cur = self.max_lengths_by_bin[chrom].get(bin_id, 0)
        self.max_lengths_by_bin[chrom][bin_id] = max(cur, feature.length)

    def get_count(self, chrom: str) -> int:
        return self.feature_counts.get(chrom, 0)

    def get_max_lengths(self, chrom: str) -> dict[int, int]:
        return self.max_lengths_by_bin.get(chrom, {})


def _callable_ref(func: Callable) -> tuple[str, str]:
    # Reject lambdas, closures, and nested defs
    if func.__name__ == "<lambda>":
        raise ValueError("feature_factory must be a top-level function or callable; lambdas/locals are not supported for pickle.")
    return (func.__module__, func.__qualname__)

def _resolve_callable(module: str, qualname: str) -> Callable:
    if "<locals>" in qualname:
        raise ValueError(f"Cannot resolve non-top-level callable: {module}.{qualname}")
    mod = importlib.import_module(module)
    obj = mod
    for part in qualname.split("."):
        obj = getattr(obj, part)
    return obj

class RegionChromosomeStore(ChromosomeFeatureStore):
    def __init__(
        self,
        chrom: str,
        feature_count: int,
        max_lengths_by_bin: dict[int, int],
        fields: Iterable[FeatureField],
        feature_factory: Callable[..., GenomicFeature],
        bin_size: int = 100_000,
    ):
        super().__init__(chromosome=chrom)
        assert feature_count > 0
        self.bin_size = bin_size
        self._starts = np.zeros(feature_count, dtype=np.uint32)
        self._ends = np.zeros(feature_count, dtype=np.uint32)
        self._feature_count = 0
        self._max_feature_length_by_bin = max_lengths_by_bin
        self._is_loaded = False

        self._fields: list[FeatureField] = list(fields)
        # Validate and record a pickle-safe reference
        self._feature_factory = feature_factory
        self._factory_ref = _callable_ref(feature_factory)

        self._extra_arrays: dict[str, np.ndarray] = {}
        for f in self._fields:
            if f.dtype == np.object_:
                self._extra_arrays[f.name] = np.empty(feature_count, dtype=np.object_)
                if f.default is not None:
                    self._extra_arrays[f.name].fill(f.default)
            else:
                self._extra_arrays[f.name] = np.zeros(feature_count, dtype=f.dtype)
                if f.default not in (None, 0):
                    self._extra_arrays[f.name].fill(f.default)

    def add(self, feature: GenomicFeature) -> None:
        if not self.index_build_mode:
            raise RuntimeError("Index build mode not active.")
        idx = self._feature_count
        self._starts[idx] = feature.start
        self._ends[idx] = feature.end
        for f in self._fields:
            self._extra_arrays[f.name][idx] = getattr(feature, f.name, f.default)
        self._feature_count += 1

    def _create_feature(self, idx: int) -> GenomicFeature:
        kwargs = {"chrom": self.chromosome, "start": int(self._starts[idx]), "end": int(self._ends[idx])}
        for f in self._fields:
            val = self._extra_arrays[f.name][idx]
            if isinstance(val, np.generic):
                val = val.item()
            kwargs[f.name] = val
        return self._feature_factory(**kwargs)

    def _binary_search_position(self, pos: int) -> int:
        left, right = 0, self._feature_count - 1
        res = self._feature_count
        while left <= right:
            mid = (left + right) // 2
            if self._starts[mid] >= pos:
                res = mid
                right = mid - 1
            else:
                left = mid + 1
        return res

    def get_by_position(self, position: int):
        if not self._is_loaded or self._feature_count == 0:
            return []
        idx = self._binary_search_position(position)
        bin_id = position // self.bin_size
        max_len = self._max_feature_length_by_bin.get(bin_id, 0)
        look_back = idx
        while look_back > 0 and position - self._starts[look_back - 1] <= max_len:
            look_back -= 1
        out = []
        for i in range(look_back, min(self._feature_count, idx + MAX_LOOK_FORWARD)):
            if self._starts[i] > position:
                break
            if self._starts[i] <= position < self._ends[i]:
                out.append(self._create_feature(i))
        return out

    def get_by_interval(self, start: int, end: int):
        if not self._is_loaded or self._feature_count == 0:
            return []
        idx = self._binary_search_position(start)
        bin_id = start // self.bin_size
        max_len = self._max_feature_length_by_bin.get(bin_id, 0)
        look_back = idx
        while look_back > 0 and start - self._starts[look_back - 1] <= max_len:
            look_back -= 1
        out = []
        for i in range(look_back, self._feature_count):
            if self._starts[i] >= end:
                break
            if self._starts[i] < end and self._ends[i] > start:
                out.append(self._create_feature(i))
        return out

    def index_build_end(self) -> None:
        super().index_build_end()
        if self._feature_count > 0:
            order = np.argsort(self._starts[:self._feature_count])
            self._starts = self._starts[order]
            self._ends = self._ends[order]
            for name, arr in self._extra_arrays.items():
                self._extra_arrays[name] = arr[order]
        self._is_loaded = True

    def __len__(self) -> int:
        return self._feature_count

    def get_features(self):
        return [self._create_feature(i) for i in range(self._feature_count)]

    # Pickle hooks
    def __getstate__(self):
        state = self.__dict__.copy()
        # Remove direct callable, keep its reference
        state["_feature_factory"] = None
        state["_factory_ref"] = self._factory_ref
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        mod, qual = self._factory_ref
        # Fail fast if someone tried to pickle a non-top-level callable
        if "<locals>" in qual:
            raise ValueError(f"feature_factory {mod}.{qual} is not top-level; define it at module scope.")
        self._feature_factory = _resolve_callable(mod, qual)
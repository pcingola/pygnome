# To do

- [x] class `GenomicFeature(BaseModel)` should be a dataclass, no a Pydantic
- [ ] GenomicFeature:
  - [ ] `MsiSiteRecord` should inherit from `GenomicFeature`
  - [ ] `Record` should inherit from `GenomicFeature` (rename to `GffRecord`)
  - [ ] `VcfRecord` should inherit from `GenomicFeature`

- [ ] `Variant` has mutiple 'ALT'???
- [ ] `VcfVariant` is useless
- [ ] `VcfRecord`: Parsing beyond chrom, start, end, re,f alt should be lazy (even spltiing tabs can take a long time for millions of samples)

- [x] `BaseSequence`: change `self._data: List[int] = []` to numpy array
- [x] Based on the concept of `DnaString`, create a new `DnaStringArray` which efficiently stores millions of dna small strings into on storage.


DnaStringArray is quite messy and doesn't work.
Let's refactor the code to create and extract string into / outof the array, then re-write the methods accordingly.
There is also code duplication between `DnaString` and `DnaStringArray` that we can avoid.

In BaseSequence, we have this which seem quite uselss:

    @property
    @abstractmethod
    def _class_name(self) -> str:
        """Name of the class for use in __repr__."""
        pass
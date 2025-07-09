# July 2025:
- Bug corrected in [extract_transp_geqdsk](lib/scripts/extract_transp_geqdsk): string and not binary was passed to the console of txpl
- [efit module](lib/python/efit/utils.py). Added compatibility with modern scipy 

# November 2024: Commit number 0688d4cab69afdb22071f8e46262612b0ac200a7
- Changed counter to int64 to avoid overflows issues
- Added gyrocenter positions to output file when `n_(p)npa > 1`
- Updated documentation

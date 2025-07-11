# July 2025 (2): b9474cbbc22a2de6e858011f2c443521af65a22d
- [Utils module](lib/python/fidasim/utils.py) write_data makes a copy of the data prior to wrtting it to avoid transposing the input data
- Added renormalization keyword to solve issue [#298](https://github.com/D3DEnergetic/FIDASIM/issues/298)
- [efit module](lib/python/efit/utils.py). Added compatibility with modern scipy 

# July 2025: 8a7207aadbd3d9bca6afc74dc6f0f5a8d50c3029
- Bug corrected in [extract_transp_geqdsk](lib/scripts/extract_transp_geqdsk): string and not binary was passed to the console of txpl
- [efit module](lib/python/efit/utils.py). Added compatibility with modern scipy 

# November 2024: Commit number 0688d4cab69afdb22071f8e46262612b0ac200a7
- Changed counter to int64 to avoid overflows issues
- Added gyrocenter positions to output file when `n_(p)npa > 1`
- Updated documentation

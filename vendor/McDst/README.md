# Bundled McDst classes

These six source files provide the ROOT classes needed by Converter and mcdst_lab:
`McEvent.h/.cxx`, `McParticle.h/.cxx`, and `McRun.h/.cxx`.

Source: [TemirB/sthbtmaker, StRoot/McDst](https://github.com/TemirB/sthbtmaker/tree/71abc2df133e0583ff9b8ee7bfd6981bfb3d011f/StRoot/McDst).
Snapshot commit: `71abc2df133e0583ff9b8ee7bfd6981bfb3d011f`.

The files are copied unchanged from that snapshot. Their ROOT class names, members,
and schema versions are preserved: McEvent 2, McParticle 2, McRun 1. This keeps the
output compatible with McDst files and readers using the same schemas.

The existing `include/ConverterMcDstLinkDef.h` generates the dictionary. No prebuilt
binaries, dictionary files, or dependency downloads are needed from sthbtmaker.

Both Make and CMake use this directory by default. To use an external McDst, set
`MCDST_SOURCE_DIR`, or `MCDST_INCLUDE_DIR` together with `MCDST_LIBRARY`.

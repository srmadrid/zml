# ZML: Zig Mathematics Library

A Zig numerical and symbolic mathematics library.

## Warning

This library is in the early stages of development and is not yet ready for use. Breaking changes are to be expected every commit, and only the most basic functionality is currently implemented.

## Current Features

- Basic N-dimensional array support.

## Installation

To use this library in your project, run

```bash
zig fetch --save git+https://github.com/srmadrid/zml
```

and add it to your build.zig file:

```zig
const zml = b.dependency("zml", .{});
exe.root_module.addImport("zml", zml.module("zml"));
```

## Notes

This library is in the early stages of development and is not yet ready for use. Breaking changes are to be expected every commit, and only the most basic functionality is currently implemented (see [Current Features](#current-features) for more information).

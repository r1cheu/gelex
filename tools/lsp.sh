#!/bin/bash
ln -s build/compile_commands.json ./


echo "CompileFlags:
  Add:
    - \"-I$(pwd)/extern/armadillo/include\"
    - \"-I$(pwd)/include\"
    " > .clangd

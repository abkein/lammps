{
  description = "C++ dev shell";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:

  let
    system = "x86_64-linux";
    pkgs = import nixpkgs { inherit system; };
  in
  {
    devShells.${system}.default = pkgs.mkShell {
      packages = with pkgs; [
        cmake
        ninja
        pkg-config
        gcc
        gdb
        gfortran
        clang

        (python3.withPackages (ps: with ps; [
          ipykernel
          pip
          bash-kernel
          ipython
          ipykernel
          jupyter
          jupyterlab
          notebook
          pyzmq
          numpy
          pandas
          scipy
          requests
          matplotlib

          adios2

          # linters
          lizard
        ]))

        openmpi
        adios2
        zlib
        blas
        lapack
        zstd

        # linters
        flawfinder
        cppcheck
      ];
    };
  };
}

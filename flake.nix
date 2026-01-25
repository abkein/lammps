{
  description = "C++ dev shell";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:

    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
      lib = pkgs.lib;
      gccToolchain = pkgs.gcc.cc; # "real" gcc, not just the wrapper
      mpi = pkgs.openmpi;
      clangTidyWrapped = pkgs.writeShellScriptBin "clang-tidy" ''
        exec ${lib.getExe' pkgs.clang "clang-tidy"} \
          --extra-arg=--gcc-toolchain=${gccToolchain} \
          --extra-arg=-stdlib=libstdc++ \
          --extra-arg=-isystem${mpi}/include \
          "$@"
      '';
      workspaceFile = pkgs.writeText "LAMMPS.code-workspace" (
        builtins.toJSON {
          folders = [
            {
              path = "/home/kein/repos/mylammps";
            }
          ];
          settings = {
            "direnv.watchForChanges" = false;
            "cmake.sourceDirectory" = "/home/kein/repos/mylammps/cmake";

            "cpplint.cpplintPath" = lib.getExe' pkgs.cpplint "cpplint";
            "cpplint.lineLength" = 300;
            "cpplint.verbose" = 0;

            "c-cpp-linter.clangTidy.path" = clangTidyWrapped;
            "c-cpp-linter.compiler.path" = lib.getExe' pkgs.clang "clang++";
            "c-cpp-linter.cppCheck.path" = lib.getExe' pkgs.cppcheck "cppcheck";
            "c-cpp-linter.general.runOnOpen" = false;
            "c-cpp-linter.general.runOnSave" = false;
            "c-cpp-linter.general.showInformationDialog" = true;
            "c-cpp-linter.general.showOutputFromLinters" = true;
          };
        }
      );
    in
    {
      devShells.${system}.default = pkgs.mkShell {
        # BETTER_CODE_VSCODE_WORKSPACE_FILE = workspaceFile;
        packages = with pkgs; [
          jq
          cmake
          ninja
          pkg-config
          gcc
          gdb
          gfortran
          clang
          # clangTidyWrapped

          (python3.withPackages (
            ps: with ps; [
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
            ]
          ))

          # libraries
          openmpi
          adios2
          fftw
          zlib
          blas
          lapack
          zstd
          gzip
          libpng

          # linters
          flawfinder
          cppcheck
          cpplint
        ];

        shellHook = ''
          export BETTER_CODE_VSCODE_WORKSPACE_FILE="$PWD/.vscode/LAMMPS.code-workspace"
          cat ${workspaceFile} | jq . >"$BETTER_CODE_VSCODE_WORKSPACE_FILE"
        '';
      };
    };
}

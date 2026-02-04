{
  description = "C++ dev shell";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:

    let
      root = "/home/kein/repos/mylammps";

      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
      lib = pkgs.lib;

      gccToolchain = pkgs.gcc.cc; # "real" gcc, not just the wrapper
      mpi = pkgs.openmpi;
      llvm = pkgs.llvmPackages_latest;

      clangTidy = lib.getExe' llvm.clang-tools "clang-tidy"; # lib.getExe' pkgs.clang "clang-tidy"
      clangTidyWrapped = pkgs.writeShellScriptBin "clang-tidy-wrapped" ''
        exec ${clangTidy} \
          --extra-arg=--gcc-toolchain=${gccToolchain} \
          --extra-arg=-stdlib=libstdc++ \
          --extra-arg=-isystem${mpi}/include \
          --extra-arg=-Qunused-arguments \
          "$@"
      '';
      vscodeDir = "${root}/.vscode";
      cppcheckBuildDir = "${root}/.cppcheck";
      cppcheckSupprPlainLoc = "${vscodeDir}/cppcheck_suppressions";
      cppcheckSupprPlain = pkgs.writeText "LAMMPS.code-workspace" ''
        noExplicitConstructor:src/nucc_cspan.hpp
      '';
      clang-tidy-conf-loc = "${vscodeDir}/.clang-tidy";
      clang-tidy-conf = pkgs.writeText "LAMMPS.clang-tidy" ''
      Checks: >
        cppcoreguidelines-pro-type-member-init

      CheckOptions:
        - key: cppcoreguidelines-pro-type-member-init.IgnoreArrays
          value: 'true'
      '';
      workspaceFileLoc = "${vscodeDir}/LAMMPS.code-workspace";
      workspaceFile = pkgs.writeText "LAMMPS.code-workspace" (
        builtins.toJSON {
          folders = [
            {
              path = root;
            }
          ];
          settings = {
            # "direnv.watchForChanges" = false;
            "cmake.sourceDirectory" = "${root}/cmake";

            "cpplint.cpplintPath" = lib.getExe' pkgs.cpplint "cpplint";
            "cpplint.lineLength" = 300;
            "cpplint.verbose" = 0;
            # "cpplint.filters" =
            #   "build/include_subdir,readability/todo,whitespace/newline,whitespace/operators,whitespace/braces,build/include_order";

            "c-cpp-linter.clangTidy.enabled" = false;
            "c-cpp-linter.clangTidy.path" = "${clangTidyWrapped}/bin/clang-tidy-wrapped";
            "c-cpp-linter.compiler.path" = lib.getExe' pkgs.clang "clang++";
            "c-cpp-linter.compiler.additionalFlags" = [
              "-Qunused-arguments"
            ];
            "c-cpp-linter.cppCheck.path" = lib.getExe' pkgs.cppcheck "cppcheck";
            "c-cpp-linter.cppCheck.additionalFlags" = [
              "--std=c++20"
              "--platform=native"
              "--check-level=exhaustive"
              "--force"
              "--cppcheck-build-dir=${cppcheckBuildDir}"
              "--inline-suppr"
              "--suppressions-list=${cppcheckSupprPlainLoc}"
              # "--enable=warning,performance,portability,information,missingInclude"
              # "--platform=unix64"
              # "-j 6"
            ];
            "c-cpp-linter.general.sourceFileExtensions" = [
              "c"
              "h"
              "cpp"
              "hpp"
            ];

            "clang-tidy.buildPath" = "${root}/build";
            "clang-tidy.lintOnSave" = false;
            "clang-tidy.executable" = clangTidy; #"${clangTidyWrapped}/bin/clang-tidy-wrapped";
            "clang-tidy.configFile" = clang-tidy-conf-loc;
            "clang-tidy.compilerArgs" = [
              "--gcc-toolchain=${gccToolchain}"
              "-stdlib=libstdc++"
              "-isystem${mpi}/include"
              "-Qunused-arguments"
            ];
            "clang-tidy.checks" = [
              "-*,boost-*,bugprone-*,concurrency-*,hicpp-*,modernize-*,performance-*,readability-*,llvm-*,misc-*,mpi-*,openmp-*"
              "-readability-magic-numbers,-readability-function-cognitive-complexity,-readability-identifier-length,-readability-math-missing-parentheses,-readability-avoid-const-params-in-decls"
              "-modernize-use-trailing-return-type,-modernize-return-braced-init-list"
              # hicpp-member-init is an alias for enabled cppcoreguidelines-pro-type-member-init
              # hicpp-special-member-functions is an alias for cppcoreguidelines-special-member-functions
              "-hicpp-signed-bitwise,-hicpp-special-member-functions,-hicpp-member-init"
              "-cppcoreguidelines-special-member-functions"
              # "-cppcoreguidelines-non-private-member-variables-in-classes"
              "-misc-non-private-member-variables-in-classes"
              "-llvm-header-guard"
              "-bugprone-easily-swappable-parameters"
            ];
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
          # clang
          llvm.clang
          llvm.clang-tools
          clangTidyWrapped

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
          mkdir -p '${vscodeDir}'
          export BETTER_CODE_VSCODE_WORKSPACE_FILE='${workspaceFileLoc}'
          cat '${workspaceFile}' | jq . > '${workspaceFileLoc}'

          mkdir -p '${cppcheckBuildDir}'
          cat '${cppcheckSupprPlain}' > '${cppcheckSupprPlainLoc}'
          cat '${clang-tidy-conf}' > '${clang-tidy-conf-loc}'
        '';
      };
    };
}

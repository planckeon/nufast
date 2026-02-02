{
  description = "nufast - Fast three-flavor neutrino oscillation probabilities in Rust";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay.url = "github:oxalica/rust-overlay";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, rust-overlay, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        rustToolchain = pkgs.rust-bin.stable.latest.default.override {
          extensions = [ "rust-src" "rust-analyzer" ];
        };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            rustToolchain
            pkg-config
            # For benchmarking
            gnuplot
            # For paper compilation
            typst
          ];

          shellHook = ''
            echo "nufast development environment"
            echo "Rust: $(rustc --version)"
            echo ""
            echo "Commands:"
            echo "  cargo build --release    # Build optimized"
            echo "  cargo bench              # Run benchmarks"
            echo "  cargo test               # Run tests"
            echo "  typst compile paper/nufast-benchmark.typ  # Compile paper"
          '';

          RUST_SRC_PATH = "${rustToolchain}/lib/rustlib/src/rust/library";
        };

        packages.default = pkgs.rustPlatform.buildRustPackage {
          pname = "nufast";
          version = "0.3.1";
          src = ./.;
          cargoLock.lockFile = ./Cargo.lock;
        };
      }
    );
}

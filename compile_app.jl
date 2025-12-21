#!/usr/bin/env julia

# Script to compile DFTatom assignment into a standalone executable
# This creates an app that can run independently of the user's Julia environment

using PackageCompiler

println("Compiling DFTatom assignment into standalone app...")
println("This may take several minutes...")

# First, ensure packages are precompiled normally
println("\nPrecompiling packages in normal environment...")
using Pkg
Pkg.precompile()

println("\nCreating standalone app...")
# Create the standalone app with incremental=true to avoid world age issues
create_app(
    ".",  # Source directory (current package)
    "build2/DFTatomApp";  # Output directory
    executables = ["dftatomcalc" => "julia_main"],  # Name => entry point
    precompile_execution_file = "src/DFTatom.jl",  # Script to precompile
    incremental = false,  # Use incremental compilation to avoid world age issues
    force = true,  # Overwrite if exists
    include_lazy_artifacts = true,  # Include all artifacts
)

println("\n" * "="^80)
println("Compilation complete!")
println("="^80)
println("\nExecutable location: build/DFTatomApp/bin/dftatomcalc")
println("\nTo run the app:")
println("  ./build/DFTatomApp/bin/dftatomcalc")
println("\nYou can copy the entire 'build/DFTatomApp' folder to any system")
println("(same OS and architecture) and run it without Julia installed!")
println("="^80)

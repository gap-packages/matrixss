Revision.init_g := "@(#)$Id$";

# Announce the package version
DeclarePackage("matrixss", "1.0", ReturnTrue);

# Install the documentation
# DeclarePackageDocumentation("matrixss", "doc");

# Read package code
ReadPkg("matrixss", "lib/code.gi");

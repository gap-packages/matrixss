Revision.("matrixss/init_g") := 
  "@(#)$Id$";

# Announce the package version
DeclarePackage("matrixss", "1.0", ReturnTrue);

# Install the documentation
# DeclarePackageDocumentation("matrixss", "doc");

# Read package declarations
ReadPkg("matrixss", "lib/code.gd");

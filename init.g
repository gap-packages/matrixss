###############################################################################
##
#W    init.g      The Matrix Schreier-Sims package                
##
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-01-10 
##
###############################################################################

# Announce the package version
#DeclarePackage("matrixss", "1.0", ReturnTrue);

# Install the documentation
#DeclarePackageDocumentation("matrixss", "doc");

# Read package declarations
ReadPkg("matrixss", "lib/code.gd");
ReadPkg("matrixss", "lib/standard.gd");
ReadPkg("matrixss", "lib/stcs.gd");
ReadPkg("matrixss", "lib/random.gd");
ReadPkg("matrixss", "lib/linear.gd");
ReadPkg("matrixss", "lib/verify.gd");
ReadPkg("matrixss", "lib/test.gd");

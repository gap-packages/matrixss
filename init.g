###############################################################################
##
#W    init.g      The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-01-10 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/init_g") := 
  "@(#)$Id$";

# Announce the package version
#DeclarePackage("matrixss", "1.0", ReturnTrue);

# Install the documentation
#DeclarePackageDocumentation("matrixss", "doc");

# Read package declarations
ReadPkg("matrixss", "lib/code.gd");
ReadPkg("matrixss", "lib/standard.gd");
ReadPkg("matrixss", "lib/stcs.gd");
ReadPkg("matrixss", "lib/random.gd");
ReadPkg("matrixss", "lib/test.gd");

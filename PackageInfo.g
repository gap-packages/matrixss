###############################################################################
##
#W    PackageInfo.g      The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-01-23
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

SetPackageInfo(rec(
        PackageName := "matrixss",
        Subtitle := "Schreier-Sims for matrix groups",
        Version := "1.0",
        #Autoload := false,
        Date := "2004-01-23",  
        PackageWWWHome := "http://matrixss.sourceforge.net/"
        ArchiveURL := 
        "http://sourceforge.net/project/showfiles.php?group_id=99072",
        ArchiveFormats := ".tar.gz,.zoo,.zip,.tar.bz2",
        Persons := 
        [rec(
             LastName := "B‰‰rnhielm",
             FirstNames := "Henrik",
             IsAuthor := true,
             IsMaintainer := true,
             Email := "redstar@linux.nu",
             WWWHome := "http://henrik.baarnhielm.net/",
             Place := "London",
             Institution := 
             "Imperial College of Science, Technology and Medicine",
             PostalAddress := 
             Concatenation(["Henrik B‰‰rnhielm\n",
                     "Imperial College, London\n",
                     "United Kingdom\n"]),
             )],
                   
#Status := "deposited",
#CommunicatedBy := "name (place)",
#AcceptDate := "MM/YYYY",
#README_URL := Concatenation( ~.PackageWWWHome, "/README" ),
#PackageInfoURL :=
#  Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
#AbstractHTML := Concatenation( [
#  "The package contains the <span class=\"pkgname\">GAP</span> ",
#  "Library of Tables of Marks"
#  ] ),
#PackageDoc := rec(
#  BookName :=
#    "TomLib",
#  ArchiveURLSubset :=
#    [ "doc", "htm" ],
#  HTMLStart :=
#    "htm/chapters.htm",
#  PDFFile :=
#    "doc/manual.pdf",
#  SixFile :=
#    "doc/manual.six",
#  LongTitle :=
#    "The GAP Library of Tables of Marks",
#  Autoload :=
#    true
#  ),
                   Dependencies := 
                    rec(
                        GAP := ">= 4.3",
                        NeededOtherPackages := [],
                        SuggestedOtherPackages := [],
                        ExternalConditions := []
                        ),
                    AvailabilityTest := ReturnTrue,
                    TestFile := "tst/test.g",
                    Keywords := ["matrix groups", "Schreier-Sims"]
                    BannerString := 
                    Concatenation(["#I loading GAP package ''matrixss'' ",
                            "in version ", ~.Version, "\n"])
                    ));

#E


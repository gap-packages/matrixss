###############################################################################
##
#W    PackageInfo.g      The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B��rnhielm
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
        Subtitle := "Schreier-Sims algorithm for matrix groups",
        Version := "0.1",
        Date := "27/07/2004",
        PackageWWWHome := "http://matrixss.sourceforge.net/",
        ArchiveURL := "http://sourceforge.net/project/showfiles.php?group_id=99072",
        PackageInfoCurrent := "http://matrixss.sourceforge.net/",
        
        ArchiveFormats := ".tar.gz,.zoo,.win.zip,.tar.bz2,.deb",
        Persons := 
        [rec(
             LastName := "B��rnhielm",
             FirstNames := "Henrik",
             IsAuthor := true,
             IsMaintainer := true,
             Email := "redstar_@sourceforge.net",
             WWWHome := "http://henrik.baarnhielm.net/",
             Place := "London",
             Institution := "Imperial College of Science, Technology and Medicine",
             PostalAddress := Concatenation(["Henrik B��rnhielm\n",
                     "Imperial College, London\n",
                     "United Kingdom\n"])
             )],
                
                Status := "alpha",
                README_URL := Concatenation( ~.PackageWWWHome, "/README" ),
                PackageInfoURL :=
                Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
                AbstractHTML := 
                Concatenation(["This package contains an implementation ",
                        "of the Schreier-Sims algorithm for matrix groups"]),
                PackageDoc := 
                rec(Archive := "",
                    BookName := "matrixss",
                    HTMLStart := "",
                    PDFFile := "",
                    SixFile := "doc/manual.six",
                    LongTitle := "",
                    Autoload := false
                    ),
                Dependencies := 
                rec(
                    GAP := ">= 4.4",
                    NeededOtherPackages := [],
                    SuggestedOtherPackages := [],
                    ExternalConditions := []
                    ),
                AvailabilityTest := ReturnTrue,
                TestFile := "tst/test.g",
                Keywords := ["matrix groups", "Schreier-Sims"],
                BannerString := 
                Concatenation(["#I loading GAP package ''matrixss'' ",
                        "in version ", ~.Version, "\n"]),
                Autoload := false,
                
                ));

###############################################################################
#E


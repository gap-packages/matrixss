###############################################################################
##
#W    PackageInfo.g      The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
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
        Version := "1.0",
        Date := "22/06/2004",
        PackageWWWHome := "http://matrixss.sourceforge.net/",
        ArchiveURL := "http://sourceforge.net/project/showfiles.php?group_id=99072",
        PackageInfoCurrent := "http://matrixss.sourceforge.net/",
        
        ArchiveFormats := ".tar.gz,.zoo,.zip,.tar.bz2",
        Persons := [rec(
                        LastName := "Bäärnhielm",
                        FirstNames := "Henrik",
                        IsAuthor := true,
                        IsMaintainer := true,
                        Email := "redstar_@sourceforge.net",
                        WWWHome := "http://henrik.baarnhielm.net/",
                        Place := "London",
                        Institution := "Imperial College of Science, Technology and Medicine",
                        PostalAddress := Concatenation(["Henrik Bäärnhielm\n",
                                "Imperial College, London\n",
                                "United Kingdom\n"])
                        )],
                   
                   Status := "dev",
                   README_URL := Concatenation( ~.PackageWWWHome, "/README" ),
                   PackageInfoURL :=
                   Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
                   AbstractHTML := 
                   Concatenation(["This package contains an implementation ",
                           "of the Schreier-Sims algorithm"]),
                   PackageDoc := rec(Archive := "",
                           BookName := "",
                           HTMLStart := "",
                           PDFFile := "",
                           SixFile := "",
                           LongTitle := "",
                           Autoload := false
                           ),
                   Dependencies := rec(
                           GAP := ">= 4.3",
                           NeededOtherPackages := [],
                           SuggestedOtherPackages := [],
                           ExternalConditions := []
                           ),
                   AvailabilityTest := ReturnTrue,
                   TestFile := "tst/test.g",
                   Keywords := ["matrix groups", "Schreier-Sims"],
                   BannerString := Concatenation(["#I loading GAP package ''matrixss'' ",
                           "in version ", ~.Version, "\n"]),
                   Autoload := false,

                   ));

#E


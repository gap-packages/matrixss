###############################################################################
##
#W    PackageInfo.g      The Matrix Schreier-Sims package
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-01-23
##

SetPackageInfo( rec(

PackageName := "matrixss",
Subtitle := "The Schreier-Sims algorithm for matrix groups",
Version := "0.9",
Date := "11/09/2004",

SourceRepository := rec(
                         Type := "git",
                         URL := Concatenation("https://github.com/gap-packages/", ~.PackageName )
                    ),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues"),
PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", ~.PackageName ),
README_URL      := Concatenation( ~.PackageWWWHome, "/README"),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                  "/releases/download/v", ~.Version,
                                  "/", ~.PackageName, "-", ~.Version ),
ArchiveFormats := ".tar.gz .tar.bz2",

Persons := [
  rec(
    LastName     := "Bäärnhielm",
    FirstNames   := "Henrik",
    IsAuthor     := true,
    IsMaintainer := true,
    Email        := "redstar_@sourceforge.net",
    WWWHome      := "http://henrik.baarnhielm.net/",
    Place        := "London",
    Institution  := Concatenation([
      "School of Mathematical Sciences, ",
      "Queen Mary, University of London",
      ])
      ),
    rec(
       IsAuthor := false,
       IsMaintainer := true,
       FirstNames := "Markus",
       LastName := "Pfeiffer",
       WWWHome := "http://www.morphism.de/~markusp/",
       Email := "markus.pfeiffer@st-andrews.ac.uk",
       PostalAddress := Concatenation( [
                                         "School of Computer Science\n",
                                         "North HaughSt Andrews\n",
                                         "Fife\n",
                                         "KY16 9SX\n",
                                         "United Kingdom" ] ),
       Place := "St Andrews",
       Institution := "University of St Andrews",
  ),
],

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
# Status := "accepted",
Status := "deposited",

AbstractHTML := "<span class=\"pkgname\">matrixss</span> is a package with \
an implementation of the Schreier-Sims algorithm for matrix groups. \
In most cases it is more efficient than the implementation in the  \
<span class=\"pkgname\">GAP</span> library.",

PackageDoc := rec(
  BookName := "matrixss",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/matrixss.html",
  PDFFile   := "doc/matrixss.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "The Schreier-Sims algorithm for matrix groups",
  Autoload := true),

Dependencies := rec(
  GAP := ">=4.4",
  NeededOtherPackages := [],
  SuggestedOtherPackages := [],
  ExternalConditions := [] ),

AvailabilityTest := ReturnTrue,
Autoload := false,
TestFile := "tst/test.g",
Keywords := ["matrix group", "Schreier-Sims"]

#BannerString := Concatenation(
#  "----------------------------------------------------------------\n",
#  "Loading matrixss ", ~.Version, "\n",
#  "by ", ~.Persons[1].FirstNames, " ", ~.Persons[1].LastName,
#        " <", ~.Persons[1].Email, ">\n",
#  "For help, type: ?matrixss package \n",
#  "----------------------------------------------------------------\n" ),
#

));

###############################################################################
#E

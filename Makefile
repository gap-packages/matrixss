# to create packages
PACKAGE_FILES = Manifest
CHECK_SUMS    = SHA1SUMS
DIR           := $(shell pwd)
CUR_DIR       := $(notdir $(DIR))

# to upload packages
FTP_USER      = anonymous
FTP_PASSWD    = redstar_@users.sourceforge.net
FTP_DIR       = /incoming
FTP_HOST      = upload.sourceforge.net

# to upload homepage
WEB_SERVER    = shell.sourceforge.net
WEB_PATH      = /home/groups/m/ma/matrixss/htdocs
SCP_USER      = redstar_
WEB_FILES     = doc/www/*.html README PackageInfo.g

# to perform test
TEST_FILE     = tst/test.g

# package default name
ifeq ($(PACKAGE),)
	PACKAGE = matrixss
endif

# output package files
PACKAGES = $(PACKAGE).tar.gz $(PACKAGE).zip $(PACKAGE).zoo $(PACKAGE).tar.bz2


all:

test:
	gap -q < $(TEST_FILE)

release: upload www

checksums:
	sha1sum --text `grep --invert-match $(CHECK_SUMS) Manifest` > $(CHECK_SUMS)

package: checksums 
	cat $(PACKAGE_FILES) | sed --expression='s/^/$(CUR_DIR)\//' | tar --verify --verbose --create --directory .. --file $(PACKAGE).tar --files-from - 
	bzip2 --keep --verbose --best $(PACKAGE).tar
	bzip2 --test --verbose $(PACKAGE).tar.bz2
	gzip --best --force --verbose $(PACKAGE).tar
	gzip --verbose --test $(PACKAGE).tar.gz
	cd .. && cat $(CUR_DIR)/$(PACKAGE_FILES) | sed --expression='s/^/$(CUR_DIR)\//' | zip -v -9 $(CUR_DIR)/$(PACKAGE).zip -@ && cd $(CUR_DIR)
	unzip -t $(PACKAGE).zip
	cd .. && cat $(CUR_DIR)/$(PACKAGE_FILES).zoo | sed --expression='s/^\([^!/]\)/$(CUR_DIR)\/\1/' | zoo achPI $(CUR_DIR)/$(PACKAGE).zoo && cd $(CUR_DIR)
	zoo xN $(PACKAGE).zoo

upload: package
	ncftpput -u $(FTP_USER) -p $(FTP_PASSWD) $(FTP_HOST) $(FTP_DIR) $(PACKAGES)

www:
	scp $(WEB_FILES) $(SCP_USER)@$(WEB_SERVER):$(WEB_PATH)

clean:
	rm --force $(PACKAGES)
	rm --recursive --force *~ *.bak

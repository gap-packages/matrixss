PACKAGE_FILES = Manifest
FTP_USER      = anonymous
FTP_PASSWD    = redstar_@users.sourceforge.net
FTP_DIR       = /incoming
FTP_HOST      = upload.sourceforge.net

ifeq ($(PACKAGE),)
	PACKAGE = matrixss
endif


package:
	tar --verify --verbose --create --directory .. --file $(PACKAGE).tar --files-from $(PACKAGE_FILES)
	gzip --best --force --verbose $(PACKAGE).tar
	gzip --verbose --test $(PACKAGE).tar.gz
	cd .. && zip -v -9 matrixss/$(PACKAGE).zip -@ < matrixss/$(PACKAGE_FILES) && cd matrixss
	unzip -t $(PACKAGE).zip
	cd .. && zoo ahPI matrixss/$(PACKAGE).zoo < matrixss/$(PACKAGE_FILES) && cd matrixss
	zoo xN $(PACKAGE).zoo

	ncftpput -u $(FTP_USER) -p $(FTP_PASSWD) $(FTP_HOST) $(FTP_DIR) $(PACKAGE).tar.gz $(PACKAGE).zip $(PACKAGE).zoo

clean:
	rm --force $(PACKAGE).tar.gz $(PACKAGE).zip $(PACKAGE).zoo
	rm --recursive --force *~ *.bak

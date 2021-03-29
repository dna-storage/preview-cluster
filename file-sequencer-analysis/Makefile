.PHONY: clean help init align cluster preview_1_set_12_16 preview_2_set_12_16 preview_3_set_12_16 preview_4_12_16

SCRIPT_DIR=scripts

init: #run init steps
	git clone https://github.com/gui11aume/starcode
	make -C starcode


clean:
	rm -rf build
	rm -rf dist
	rm -rf dnastorage.egg-info
	rm -rf generate.egg-info

install:
	python setup.py install --user


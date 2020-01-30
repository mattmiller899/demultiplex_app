SINGULARITY = /usr/local/bin/singularity
IMG_NAME = demultiplexer.img
DEF_NAME = image.def

run:
	sudo $(SINGULARITY) run $(IMG_NAME)

clean:
	rm -f $(IMG_NAME)

img: clean
	sudo $(SINGULARITY) build $(IMG_NAME) $(DEF_NAME)

shell:
	sudo $(SINGULARITY) shell --writable -B $(shell pwd):/home $(IMG_NAME)

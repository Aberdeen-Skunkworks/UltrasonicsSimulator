all:
	$(MAKE) main

main:
	$(MAKE) -C build/

clean:
	if [ -d "build/" ]; then \
		rm -r build; \
	fi
	mkdir build
	cd build; \
	cmake ..
	$(MAKE) main

simple_commit:
	git commit -a -m"more work"

simple_push:
	make simple_commit
	git push

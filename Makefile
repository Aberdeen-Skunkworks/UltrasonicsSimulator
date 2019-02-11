all:
	$(MAKE) main

main:
	$(MAKE) -C build/

clean:
	rm -r build
	mkdir build
	cd build; \
	cmake ..
	$(MAKE) main

simple_commit:
	git commit -a -m"more work"


simple_push:
	make simple_commit
	git push

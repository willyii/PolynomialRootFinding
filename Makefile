.PHONY: build
build:
	mkdir -p build
	cd build && \
	cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 .. && \
	mv compile_commands.json ../src/ && \
	make

.PHONY: index
index:
	mkdir -p build
	cd build && \
	cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 .. && \
	mv compile_commands.json ../src/

.PHONY: debug
debug:
	mkdir -p build
	cd build && \
	cmake -DCMAKE_BUILD_TYPE=debug .. && \
	make

.PHONY: clean
clean:
	rm -rf build

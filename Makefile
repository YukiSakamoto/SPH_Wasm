all:
	emcc em.cpp -o em.js -s WASM=1 -O3 \
		-sMODULARIZE -sEXPORT_NAME=createModule  \
		--bind -s ENVIRONMENT=web \
		-s EXPORT_ES6=1 \
		-s ALLOW_MEMORY_GROWTH=1 \
		-s EXPORTED_RUNTIME_METHODS=['HEAPF32'] \
		-msimd128
debug:
	emcc em.cpp -o em.js -s WASM=1 -O0 \
		-sMODULARIZE -sEXPORT_NAME=createModule  \
		--bind -s ENVIRONMENT=web \
		-s EXPORT_ES6=1 \
		-s ALLOW_MEMORY_GROWTH=1 \
		-s EXPORTED_RUNTIME_METHODS=['HEAPF32'] \
		-msimd128 \
		-g3

out:
	emcc em.cpp -o em.o -s WASM=1 -O3 \
		-sMODULARIZE -sEXPORT_NAME=createModule  \
		--bind -s ENVIRONMENT=web \
		-s EXPORT_ES6=1 \
		-s ALLOW_MEMORY_GROWTH=1 \
		-s EXPORTED_RUNTIME_METHODS=['HEAPF32'] \
		-msimd128

native:
	g++ em.cpp -o em.x  -O3 
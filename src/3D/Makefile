export SERVER = 1
export OPT=-march=native
include ./rbnd/src/Makefile.defs

CC = mpicc

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
.PHONY: all librebound
all: 3D.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM)
	@echo ""
	@echo "Compilation successful. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@echo ""

librebound:
	@echo "Compiling shared library $(LIBREBOUND) ..."
	$(MAKE) -C ./rbnd/src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)
	@echo ""

clean:
	@echo "Cleaning up shared library $(LIBREBOUND) ..."
	$(MAKE) -C ./rbnd/src/ clean
	@echo "Cleaning up local directory ..."
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)

rebound_webgl.html: 3D.c
	@echo "Compiling 3D.c with emscripten (WebGL enabled)..."
	emcc -O3 -I./rbnd/src/ ./rbnd/src/*.c 3D.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ./rbnd/web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: 3D.c
	@echo "Compiling 3D.c with emscripten (WebGL disabled)..."
	emcc -O3 -I./rbnd/src/ ./rbnd/src/*.c 3D.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ./rbnd/web_client/shell_rebound_console.html -o rebound_console.html
